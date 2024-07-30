library(tidyverse)
library(rstan)

options(mc.cores = parallel::detectCores())

set.seed(604)

# Again, note that this code is just reproducing https://khakieconomics.github.io/2018/12/27/Ranked-random-coefficients-logit.html,
# although with a more real-world feeling number of respondents and
# and total choices. If this is taking uncomfortably long, feel free to scale
# down I, the results should reproduce just fine.

# Number of individuals
I <- 100
# Number of tasks per individual
Tasks <- 10
# Number of choices per task
J <- 5
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

# Create X -- let's make this a dummy matrix
X <- matrix(sample(0:1, Tasks*I*J*P, replace = T), Tasks*I*J, P)
# Each of the rows in this matrix correspond to a choice presented to a given individual
# in a given task

indexes <- crossing(individual = 1:I, task = 1:Tasks, option = 1:J) %>% 
  mutate(row = 1:n())

# Write a Gumbel random number generator using inverse CDF trick
rgumbel <- function(n, mu = 0, beta = 1) mu - beta * log(-log(runif(n)))
mean(rgumbel(1e6))

# Ok, now we need to simulate choices. Each person in each task compares each 
# choice according to X*beta_i + epsilon, where epsilon is gumbel distributed. 
# They return their rankings. 

# Modify the ranked_options data frame to include all necessary information
ranked_options <- indexes %>% 
  group_by(individual, task) %>% 
  mutate(
    fixed_utility = as.numeric(X[row,] %*% as.numeric(beta_i[first(individual),])),
    plus_gumbel_error = fixed_utility + rgumbel(n()),
    true_rank = rank(-plus_gumbel_error),
    true_order = order(true_rank),
    observed_order = case_when(
      true_order == 1 ~ J,  # Worst choice
      true_order == J ~ 1,  # Best choice
      TRUE ~ 3  # Tie all middle orders
    ),
    best_choice = as.numeric(true_order == J),
    worst_choice = as.numeric(true_order == 1)
  )

tt <- ranked_options %>% 
  group_by(individual, task) %>%
  summarise(start = min(row), 
            end = max(row)) %>% 
  ungroup %>%
  mutate(task_number = 1:n())

stan_data <- list(
  N = nrow(X),
  T = nrow(tt),
  I = I, 
  P = P, 
  P2 = P2, 
  K = J, 
  rank_order = ranked_options$observed_order,
  X = X, 
  X2 = W, 
  task = tt$task_number, 
  task_individual = tt$individual,
  start = tt$start, 
  end = tt$end
)


efron_simplified_rol <- "
functions {
  real rank_logit_ties_lpmf(int[] rank_order, vector delta) {
    int K = rows(delta);
    vector[K] tmp = delta[rank_order];
    real out = 0;
    int current_pos = 1;

    while (current_pos <= K) {
      int tied_count = 1;
      
      // Find ties
      while (current_pos + tied_count <= K && rank_order[current_pos] == rank_order[current_pos + tied_count]) {
        tied_count += 1;
      }

      if (tied_count == 1) {
        // No tie, use original calculation
        out += tmp[current_pos] - log_sum_exp(tmp[current_pos:]);
      } else {
        // Tied (unknown) orderings
        out += log_sum_exp(tmp[current_pos:(current_pos+tied_count-1)]) - log_sum_exp(tmp[current_pos:]);
      }

      current_pos += tied_count;
    }
    
    return out;
  }
}

data {
  int N; // number of rows
  int T; // number of individual-choice sets/task combinations
  int I; // number of Individuals
  int P; // number of covariates that vary by choice
  int P2; // number of covariates that vary by individual
  int K; // number of choices
  
  int rank_order[N]; // The vector describing the index (within each task) of the first, second, third, ... choices. 
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
  matrix[I, P] beta_individual = rep_matrix(beta', I) + X2 * Gamma' + z * diag_pre_multiply(tau, L_Omega);
}

model {
  tau ~ normal(0, .5);
  beta ~ normal(0, 1);
  to_vector(z) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(Gamma) ~ normal(0, 1);

  for(t in 1:T) {
    vector[K] utilities;
    utilities = X[start[t]:end[t]] * beta_individual[task_individual[t]]';
    target += rank_logit_ties_lpmf(rank_order[start[t]:end[t]] | utilities);
  }
}
"

# Compile the model
compiled_model <- stan_model(model_code = efron_simplified_rol)

# Fit the model
fit <- sampling(compiled_model, 
                data = stan_data, 
                chains = 4, 
                iter = 2000, 
                warmup = 1000,
                cores = 4)

summary(fit)

rol_choice <- as.data.frame(fit, pars = "beta_individual") %>%
  gather(Parameter, Value) %>%
  group_by(Parameter) %>%
  summarise(median = median(Value),
            lower = quantile(Value, .05),
            upper = quantile(Value, .95)) %>%
  mutate(individual = str_extract(Parameter, "[0-9]+(?=,)") %>% parse_number,
         column = str_extract(Parameter, ",[0-9]{1,2}") %>% parse_number) %>%
  arrange(individual, column) %>%
  mutate(`True value` = as.numeric(t(beta_i)),
         Dataset = "Rank-Ordered Logit w ties")

ggplot(rol_choice, aes(x = `True value`, y = median, color = Dataset)) +
  geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_point(aes(y = median), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(y = "Utility Estimates",
       title = "Utility Estimates from the Three Models",
       subtitle = "With interior 90% credibility intervals") +
  scale_color_manual(values = c("Rank-Ordered Logit" = "green")) +
  facet_wrap(~Dataset) +
  theme(legend.position="bottom")


# Debugging

# try ROL to see if we learn anything
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

data_list_ranked_rcl <- list(N = nrow(X),
                             T = nrow(tt),
                             I = I, 
                             P = P, 
                             P2 = P2, 
                             K = J, 
                             # NOTE!! This is the tricky bit -- we use the order of the ranks (within task)
                             # Not the raw rank orderings. This is how we get the likelihood evaluation to be pretty quick
                             rank_order = ranked_options$true_order,
                             X = X, 
                             X2 = W, 
                             task = tt$task_number, 
                             task_individual = tt$individual,
                             start = tt$start, 
                             end = tt$end)

# Compile the model
rol_model <- stan_model(model_code = ranked)

# Fit the model
fit_rol <- sampling(rol_model, data = data_list_ranked_rcl, 
                    iter = 2000, warmup = 1000, chains = 4, cores = 4)

# Check convergence
rol_choice <- as.data.frame(fit_rol, pars = "beta_individual") %>%
  gather(Parameter, Value) %>%
  group_by(Parameter) %>%
  summarise(median = median(Value),
            lower = quantile(Value, .05),
            upper = quantile(Value, .95)) %>%
  mutate(individual = str_extract(Parameter, "[0-9]+(?=,)") %>% parse_number,
         column = str_extract(Parameter, ",[0-9]{1,2}") %>% parse_number) %>%
  arrange(individual, column) %>%
  mutate(`True value` = as.numeric(t(beta_i)),
         Dataset = "Rank-Ordered Logit w ties")

ggplot(rol_choice, aes(x = `True value`, y = median, color = Dataset)) +
  geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_point(aes(y = median), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(y = "Utility Estimates",
       title = "Utility Estimates from the Three Models",
       subtitle = "With interior 90% credibility intervals") +
  scale_color_manual(values = c("Rank-Ordered Logit" = "green")) +
  facet_wrap(~Dataset) +
  theme(legend.position="bottom")
