library(tidyverse)
library(rstan)
library(ibd)
options(mc.cores = parallel::detectCores())

set.seed(604)

# Again, note that this code is just reproducing https://khakieconomics.github.io/2018/12/27/Ranked-random-coefficients-logit.html,
# although with a more real-world feeling number of respondents and
# and total choices. If this is taking uncomfortably long, feel free to scale
# down I, the results should reproduce just fine.

# Number of individuals
I <- 30
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

bibd_result <- bibd(v = 10, b = 120, r = 36, k = 3, lambda = 9)

choice_sets <- t(bibd_result$design)


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


compiled_best_choice_model <- stan_model(model_code = best)

best_choice_fit <- sampling(compiled_best_choice_model, 
                            data = data_list_best_choice, 
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


