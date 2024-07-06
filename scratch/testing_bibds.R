library(tidyverse)
library(rstan)
library(ibd)
library(crossdes)
options(mc.cores = parallel::detectCores())

set.seed(604)

# Again, note that this code is just reproducing https://khakieconomics.github.io/2018/12/27/Ranked-random-coefficients-logit.html,
# although with a more real-world feeling number of respondents and
# and total choices. If this is taking uncomfortably long, feel free to scale
# down I, the results should reproduce just fine.

# Number of individuals
I <- 800
# Number of tasks per individual
Tasks <- 10
# Number of choices per task
J <- 3
# Dimension of covariate matrix
P <- 5
# Dimension of demographic matrix
P2 <- 6

# demographic matrix
W <- matrix(rnorm(I*P2), I, P2)
# Loading matrix
Gamma <- matrix(rnorm(P*P2), P, P2)

# Show W * t(Gamma) to make sure it looks right
W %*% t(Gamma)

# Correlation of decisionmaker random slopes
Omega <- cor(matrix(rnorm(P*(P+2)), P+2, P))

# Scale of decisionmaker random slopes
tau <- abs(rnorm(P, 0, .5))

# Covariance matrix of decisionmaker random slopes
Sigma <- diag(tau) %*% Omega %*% diag(tau)

# Centers of random slopes
beta <- rnorm(P)

# Individual slopes
beta_i <- MASS::mvrnorm(I, beta, Sigma) + W %*% t(Gamma)

# Again, quick plot to sanity check
plot(as.data.frame(beta_i))

bibd <- find.BIB(trt = Tasks, b = I*Tasks, k = J)

# Checks if the design is balanced wrt to both rows and columns.
# The above function is not guaranteed to produce a valid BIBD; it may produce
# a close but imperfect design, which is sufficient for our purposes. I will
# however confirm it's connected, since we'd start to lose quite a lot of efficiency
# if it were not. 

# Returns 1 if the design is connected
is.connected(bibd)


# Create a mapping from BIBD options to attribute levels; idea is to make
# identically shaped inputs to our previous test, but using the conditions
# from the BIBD instead of sampling randomly for X.
option_attributes <- matrix(sample(0:1, 10*P, replace = TRUE), nrow = 10, ncol = P)

# Function to convert BIBD task to attribute matrix
bibd_to_attributes <- function(bibd_row) {
  matrix(c(
    option_attributes[bibd_row[1], ],
    option_attributes[bibd_row[2], ],
    option_attributes[bibd_row[3], ]
  ), nrow = 3, byrow = TRUE)
}

# Generate X matrix based on BIBD
X <- do.call(rbind, lapply(1:nrow(bibd), function(i) bibd_to_attributes(bibd[i,])))

indexes <- crossing(individual = 1:I, task = 1:Tasks, option = 1:J) %>% 
  mutate(row = 1:n())

# Write a Gumbel random number generator using inverse CDF trick
rgumbel <- function(n, mu = 0, beta = 1) mu - beta * log(-log(runif(n)))
mean(rgumbel(1e6))

# Ok, now we need to simulate choices. Each person in each task compares each 
# choice according to X*beta_i + epsilon, where epsilon is gumbel distributed. 
# They return their rankings. 

ranked_options <- indexes %>% 
  group_by(individual, task) %>% 
  mutate(fixed_utility = as.numeric(X[row,] %*% as.numeric(beta_i[first(individual),])),
         plus_gumbel_error = fixed_utility + rgumbel(n()),
         rank = rank(-plus_gumbel_error),
         # We're going to use the order rather than the rank in the Stan part of the model
         order = order(rank),
         # And here we create a dummy vector for the best choice
         best_choice = as.numeric(1:n() == which.max(plus_gumbel_error)),
         worst_choice = as.numeric(1:n() == which.min(plus_gumbel_error))
  )


tt <- ranked_options %>% 
  group_by(individual, task) %>%
  summarise(start = min(row), 
            end = max(row)) %>% 
  ungroup %>%
  mutate(task_number = 1:n())

best <- "// saved as mixed_conditional_individual_effects.stan
data {
  int N; // number of rows
  int T; // number of inidvidual-choice sets/task combinations
  int I; // number of Individuals
  int P; // number of covariates that vary by choice
  int P2; // number of covariates that vary by individual
  int K; // number of choices
  
  vector<lower = 0, upper = 1>[N] choice; // binary indicator for choice
  matrix[N, P] X; // choice attributes
  matrix[I, P2] X2; // individual attributes
  
  int task[T]; // index for tasks
  int task_individual[T]; // index for individual
  int start[T]; // the starting observation for each task
  int end[T]; // the ending observation for each task
}
parameters {
  vector[P] beta; // hypermeans of the part-worths
  matrix[P, P2] Gamma; // coefficient matrix on individual attributes
  vector<lower = 0>[P] tau; // diagonal of the part-worth covariance matrix
  matrix[I, P] z; // individual random effects (unscaled)
  cholesky_factor_corr[P] L_Omega; // the cholesky factor of the correlation matrix of tastes/part-worths
}
transformed parameters {
  // here we use the reparameterization discussed on slide 30
  matrix[I, P] beta_individual = rep_matrix(beta', I) + X2 * Gamma' + z*diag_pre_multiply(tau, L_Omega);
}
model {
  // create a temporary holding vector
  vector[N] log_prob;
  
  // priors on the parameters
  tau ~ normal(0, .5);
  beta ~ normal(0, .5);
  to_vector(z) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(Gamma) ~ normal(0, 1);
  
  // log probabilities of each choice in the dataset
  for(t in 1:T) {
    vector[K] utilities; // tmp vector holding the utilities for the task/individual combination
    // add utility from product attributes with individual part-worths/marginal utilities
    utilities = X[start[t]:end[t]]*beta_individual[task_individual[t]]';
    
    log_prob[start[t]:end[t]] = log_softmax(utilities);
  }
  
  // use the likelihood derivation on slide 29
  target += log_prob' * choice;
}"

# MaxDiff model
best_worst <- "// saved as mixed_conditional_individual_effects.stan
data {
  int N; // number of rows
  int T; // number of inidvidual-choice sets/task combinations
  int I; // number of Individuals
  int P; // number of covariates that vary by choice
  int P2; // number of covariates that vary by individual
  int K; // number of choices
  
  vector<lower = 0, upper = 1>[N] choice; // binary indicator for choice
  vector<lower = 0, upper = 1>[N] worst_choice; // binary indicator for worst choice
  matrix[N, P] X; // choice attributes
  matrix[I, P2] X2; // individual attributes
  
  int task[T]; // index for tasks
  int task_individual[T]; // index for individual
  int start[T]; // the starting observation for each task
  int end[T]; // the ending observation for each task
}
parameters {
  vector[P] beta; // hypermeans of the part-worths
  matrix[P, P2] Gamma; // coefficient matrix on individual attributes
  vector<lower = 0>[P] tau; // diagonal of the part-worth covariance matrix
  matrix[I, P] z; // individual random effects (unscaled)
  cholesky_factor_corr[P] L_Omega; // the cholesky factor of the correlation matrix of tastes/part-worths
}
transformed parameters {
  // here we use the reparameterization discussed on slide 30
  matrix[I, P] beta_individual = rep_matrix(beta', I) + X2 * Gamma' + z*diag_pre_multiply(tau, L_Omega);
}
model {
  // create a temporary holding vector
  vector[N] log_prob;
  vector[N] log_prob_worst;
  
  // priors on the parameters
  tau ~ normal(0, .5);
  beta ~ normal(0, .5);
  to_vector(z) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(Gamma) ~ normal(0, 1);
  
  // log probabilities of each choice in the dataset
  for(t in 1:T) {
    vector[K] utilities; // tmp vector holding the utilities for the task/individual combination
    // add utility from product attributes with individual part-worths/marginal utilities
    utilities = X[start[t]:end[t]]*beta_individual[task_individual[t]]';
    
    log_prob[start[t]:end[t]] = log_softmax(utilities);
    log_prob_worst[start[t]:end[t]] = log_softmax(-utilities);
  }
  
  // use the likelihood derivation on slide 29
  target += log_prob' * choice;
  target += log_prob_worst' * worst_choice;
}"

# Rank-Ordered Logit model
ranked <- "// saved as ranked_rcl.stan
functions {
  real rank_logit_lpmf(int[] rank_order, vector delta) {
    // We reorder the raw utilities so that the first rank is first, second rank second... 
    vector[rows(delta)] tmp = delta[rank_order];
    real out;
    // ... and sequentially take the log of the first element of the softmax applied to the remaining
    // unranked elements.
    for(i in 1:(rows(tmp) - 1)) {
      if(i == 1) {
        out = tmp[1] - log_sum_exp(tmp);
      } else {
        out += tmp[i] - log_sum_exp(tmp[i:]);
      }
    }
    // And return the log likelihood of observing that ranking
    return(out);
  }
}
data {
  int N; // number of rows
  int T; // number of inidvidual-choice sets/task combinations
  int I; // number of Individuals
  int P; // number of covariates that vary by choice
  int P2; // number of covariates that vary by individual
  int K; // number of choices
  
  int rank_order[N]; // The vector describing the index (within each task) of the first, second, third, ... choices. 
  // In R, this is order(-utility) within each task
  matrix[N, P] X; // choice attributes
  matrix[I, P2] X2; // individual attributes
  
  int task[T]; // index for tasks
  int task_individual[T]; // index for individual
  int start[T]; // the starting observation for each task
  int end[T]; // the ending observation for each task
}
parameters {
  vector[P] beta; // hypermeans of the part-worths
  matrix[P, P2] Gamma; // coefficient matrix on individual attributes
  vector<lower = 0>[P] tau; // diagonal of the part-worth covariance matrix
  matrix[I, P] z; // individual random effects (unscaled)
  cholesky_factor_corr[P] L_Omega; // the cholesky factor of the correlation matrix of tastes/part-worths
}
transformed parameters {
  // here we use the reparameterization discussed on slide 30
  matrix[I, P] beta_individual = rep_matrix(beta', I) + X2 * Gamma' + z * diag_pre_multiply(tau, L_Omega);
}
model {
  // priors on the parameters
  tau ~ normal(0, .5);
  beta ~ normal(0, 1);
  to_vector(z) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(Gamma) ~ normal(0, 1);
  
  // log probabilities of each choice in the dataset
  for(t in 1:T) {
    vector[K] utilities; // tmp vector holding the utilities for the task/individual combination
    // add utility from product attributes with individual part-worths/marginal utilities
    utilities = X[start[t]:end[t]]*beta_individual[task_individual[t]]';
    rank_order[start[t]:end[t]] ~ rank_logit(utilities);
  }
}"

data_list_best_choice <-list(N = nrow(X),
                             T = nrow(tt),
                             I = I, 
                             P = P, 
                             P2 = P2, 
                             K = J, 
                             # NOTE!! This is the tricky bit -- we use the order of the ranks (within task)
                             # Not the raw rank orderings. This is how we get the likelihood evaluation to be pretty quick
                             choice = ranked_options$best_choice,
                             X = X, 
                             X2 = W, 
                             task = tt$task_number, 
                             task_individual = tt$individual,
                             start = tt$start, 
                             end = tt$end)

data_list_best_worst <-list(N = nrow(X),
                            T = nrow(tt),
                            I = I, 
                            P = P, 
                            P2 = P2, 
                            K = J,
                            choice = ranked_options$best_choice,
                            worst_choice = ranked_options$worst_choice,
                            X = X, 
                            X2 = W, 
                            task = tt$task_number, 
                            task_individual = tt$individual,
                            start = tt$start, 
                            end = tt$end)

data_list_ranked_rcl <- list(N = nrow(X),
                             T = nrow(tt),
                             I = I, 
                             P = P, 
                             P2 = P2, 
                             K = J, 
                             # NOTE!! This is the tricky bit -- we use the order of the ranks (within task)
                             # Not the raw rank orderings. This is how we get the likelihood evaluation to be pretty quick
                             rank_order = ranked_options$order,
                             X = X, 
                             X2 = W, 
                             task = tt$task_number, 
                             task_individual = tt$individual,
                             start = tt$start, 
                             end = tt$end)

compiled_best_choice_model <- stan_model(model_code = best)
best_choice_fit <- sampling(compiled_best_choice_model, 
                            data = data_list_best_choice, 
                            iter = 800)

# Compile and fit MaxDiff model
compiled_maxdiff_model <- stan_model(model_code = best_worst)
maxdiff_fit <- sampling(compiled_maxdiff_model,
                        data = data_list_best_worst,
                        iter = 800)

# Compile and fit ROL model
compiled_rol_model <- stan_model(model_code = ranked)
rol_fit <- sampling(compiled_rol_model,
                    data = data_list_ranked_rcl,
                    iter = 800)

best_choice <- as.data.frame(best_choice_fit, pars = "beta_individual") %>%
  gather(Parameter, Value) %>%
  group_by(Parameter) %>%
  summarise(median = median(Value),
            lower = quantile(Value, .05),
            upper = quantile(Value, .95)) %>%
  mutate(individual = str_extract(Parameter, "[0-9]+(?=,)") %>% parse_number,
         column = str_extract(Parameter, ",[0-9]{1,2}") %>% parse_number) %>%
  arrange(individual, column) %>%
  mutate(`True value` = as.numeric(t(beta_i)),
         Dataset = "Best Choice")

# Process MaxDiff results
maxdiff_choice <- as.data.frame(maxdiff_fit, pars = "beta_individual") %>%
  gather(Parameter, Value) %>%
  group_by(Parameter) %>%
  summarise(median = median(Value),
            lower = quantile(Value, .05),
            upper = quantile(Value, .95)) %>%
  mutate(individual = str_extract(Parameter, "[0-9]+(?=,)") %>% parse_number,
         column = str_extract(Parameter, ",[0-9]{1,2}") %>% parse_number) %>%
  arrange(individual, column) %>%
  mutate(`True value` = as.numeric(t(beta_i)),
         Dataset = "MaxDiff")

# Process ROL results
rol_choice <- as.data.frame(rol_fit, pars = "beta_individual") %>%
  gather(Parameter, Value) %>%
  group_by(Parameter) %>%
  summarise(median = median(Value),
            lower = quantile(Value, .05),
            upper = quantile(Value, .95)) %>%
  mutate(individual = str_extract(Parameter, "[0-9]+(?=,)") %>% parse_number,
         column = str_extract(Parameter, ",[0-9]{1,2}") %>% parse_number) %>%
  arrange(individual, column) %>%
  mutate(`True value` = as.numeric(t(beta_i)),
         Dataset = "Rank-Ordered Logit")

ggplot(best_choice, aes(x = `True value`, y = median, color = Dataset)) +
  geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_point(aes(y = median), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(y = "Utility Estimates",
       title = "Utility Estimates from the Two Models",
       subtitle = "With interior 90% credibility intervals") +
  scale_color_manual(values = c("Best Choice" = "blue")) +
  facet_wrap(~Dataset) +
  theme(legend.position="none")

# Combine all results
all_results <- bind_rows(best_choice, maxdiff_choice, rol_choice)

# Plot combined results
ggplot(all_results, aes(x = `True value`, y = median, color = Dataset)) +
  geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_point(aes(y = median), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(y = "Utility Estimates",
       title = "Utility Estimates from the Three Models",
       subtitle = "With interior 90% credibility intervals") +
  scale_color_manual(values = c("Best Choice" = "blue", "MaxDiff" = "red", "Rank-Ordered Logit" = "green")) +
  facet_wrap(~Dataset) +
  theme(legend.position="bottom")

all_results %>%
  group_by(Dataset) %>%
  summarize(RMSE = sqrt(mean((`True value` - median)^2))
