library(tidyverse)
library(rstan)
library(combinat)

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
      true_order == J ~ J,  # Worst choice
      true_order == 1 ~ 1,  # Best choice
      TRUE ~ 3  # Tie all middle orders
    ),
    best_choice = as.numeric(true_order == 1),
    worst_choice = as.numeric(true_order == J)
  )

tt <- ranked_options %>% 
  group_by(individual, task) %>%
  summarise(start = min(row), 
            end = max(row)) %>% 
  ungroup %>%
  mutate(task_number = 1:n())

n_tied <- J - 2  # Adjust this to match your DGP
permutations <- permn(n_tied)
permutation_matrix <- do.call(rbind, permutations)
n_permutations <- nrow(permutation_matrix)

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
  end = tt$end,
  n_tied = n_tied,
  permutations = permutation_matrix,
  n_permutations = n_permutations
)


efron_simplified_rol <- "
functions {
  real rank_logit_ties_lpmf(int[] y, vector delta, int[,] permutations, int n_permutations, int n_tied) {
    int K = rows(delta);
    vector[K] sorted_delta = delta[y];
    real out = 0;
    
    // Handle known best
    out += sorted_delta[1] - log_sum_exp(sorted_delta);
    
    // Compute normalizing factor once
    real normalizing_factor = log_sum_exp(sorted_delta[2:K]);
    
    // Handle tied middle ranks
    real perm_sum = 0;
    for (p in 1:n_permutations) {
      real perm_ll = 0;
      for (i in 1:n_tied) {
        int idx = permutations[p, i];
        perm_ll += sorted_delta[1+idx];
      }
      perm_sum += exp(perm_ll - n_tied * normalizing_factor);
    }
    out += log(perm_sum / n_permutations);
    
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
  
  int<lower=2> n_tied;
  int<lower=1> n_permutations;
  int permutations[n_permutations, n_tied];
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
    rank_order[start[t]:end[t]] ~ rank_logit_ties(utilities, permutations, n_permutations, n_tied);
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

normalize_utilities <- function(utilities) {
  (utilities - mean(utilities)) / sd(utilities)
}

rol_choice <- as.data.frame(fit, pars = "beta_individual") %>%
  gather(Parameter, Value) %>%
  mutate(
    individual = str_extract(Parameter, "[0-9]+(?=,)") %>% parse_number(),
    column = str_extract(Parameter, ",[0-9]{1,2}") %>% parse_number()
  ) %>%
  group_by(individual) %>%
  mutate(Value_normalized = normalize_utilities(Value)) %>%
  group_by(individual, column) %>%
  summarise(
    median = median(Value_normalized),
    lower = quantile(Value_normalized, 0.05),
    upper = quantile(Value_normalized, 0.95),
    .groups = 'drop'
  ) %>%
  ungroup() %>%
  mutate(
    True_value = as.numeric(t(beta_i)),
    Dataset = "Rank-Ordered Logit w ties"
  ) %>%
  group_by(individual) %>%
  mutate(True_value_normalized = normalize_utilities(True_value)) %>%
  ungroup()

ggplot(rol_choice, aes(x = True_value_normalized, y = median, color = Dataset)) +
  geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_point(aes(y = median), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Normalized True Utility",
       y = "Normalized Estimated Utility",
       title = "Normalized Utility Estimates vs True Values",
       subtitle = "With interior 90% credibility intervals") +
  scale_color_manual(values = c("Rank-Ordered Logit w ties" = "green")) +
  theme_minimal() +
  coord_fixed(ratio = 1)

rmse <- sqrt(mean((rol_choice$True_value_normalized - rol_choice$median)^2))
print(paste("RMSE:", round(rmse, 4)))

# Compare to Best-worst
stan_data_maxdiff <- list(
  N = nrow(X),
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
  end = tt$end
)

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

compiled_maxdiff_model <- stan_model(model_code = best_worst)

maxdiff_fit <- sampling(compiled_maxdiff_model, 
                        data = stan_data_maxdiff, 
                        chains = 4, 
                        iter = 800, 
                        warmup = 400,
                        cores = 4)

# Process MaxDiff results
maxdiff_choice <- as.data.frame(maxdiff_fit, pars = "beta_individual") %>%
  gather(Parameter, Value) %>%
  mutate(
    individual = str_extract(Parameter, "[0-9]+(?=,)") %>% parse_number(),
    column = str_extract(Parameter, ",[0-9]{1,2}") %>% parse_number()
  ) %>%
  group_by(individual) %>%
  mutate(Value_normalized = normalize_utilities(Value)) %>%
  group_by(individual, column) %>%
  summarise(
    median = median(Value_normalized),
    lower = quantile(Value_normalized, 0.05),
    upper = quantile(Value_normalized, 0.95),
    .groups = 'drop'
  ) %>%
  ungroup() %>%
  mutate(
    True_value = as.numeric(t(beta_i)),
    Dataset = "MaxDiff"
  ) %>%
  group_by(individual) %>%
  mutate(True_value_normalized = normalize_utilities(True_value)) %>%
  ungroup()

# Combine results
combined_results <- bind_rows(rol_choice, maxdiff_choice)

# Create comparison plot
ggplot(combined_results, aes(x = True_value_normalized, y = median, color = Dataset)) +
  geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_point(aes(y = median), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Normalized True Utility",
       y = "Normalized Estimated Utility",
       title = "Normalized Utility Estimates vs True Values",
       subtitle = "With interior 90% credibility intervals") +
  scale_color_manual(values = c("Rank-Ordered Logit w ties" = "green", "MaxDiff" = "blue")) +
  theme_minimal() +
  coord_fixed(ratio = 1) +
  facet_wrap(~Dataset)

# Calculate RMSE for both models
rmse_rol <- sqrt(mean((rol_choice$True_value_normalized - rol_choice$median)^2))
rmse_maxdiff <- sqrt(mean((maxdiff_choice$True_value_normalized - maxdiff_choice$median)^2))

# Print RMSE values
print(paste("RMSE Rank-Ordered Logit with ties:", round(rmse_rol, 4)))
print(paste("RMSE MaxDiff:", round(rmse_maxdiff, 4)))

# Add RMSE to the plot
plot_with_rmse <- last_plot() +
  labs(subtitle = paste("With interior 90% credibility intervals\n",
                        "RMSE ROL:", round(rmse_rol, 4), 
                        "RMSE MaxDiff:", round(rmse_maxdiff, 4)))

print(plot_with_rmse)

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

summary(fit_rol)

# Check convergence
rol_choice_no_ties <- as.data.frame(fit_rol, pars = "beta_individual") %>%
  gather(Parameter, Value) %>%
  mutate(
    individual = str_extract(Parameter, "[0-9]+(?=,)") %>% parse_number(),
    column = str_extract(Parameter, ",[0-9]{1,2}") %>% parse_number()
  ) %>%
  group_by(individual) %>%
  mutate(Value_normalized = normalize_utilities(Value)) %>%
  group_by(individual, column) %>%
  summarise(
    median = median(Value_normalized),
    lower = quantile(Value_normalized, 0.05),
    upper = quantile(Value_normalized, 0.95),
    .groups = 'drop'
  ) %>%
  ungroup() %>%
  mutate(
    True_value = as.numeric(t(beta_i)),
    Dataset = "Rank-Ordered Logit (no ties)"
  ) %>%
  group_by(individual) %>%
  mutate(True_value_normalized = normalize_utilities(True_value)) %>%
  ungroup()

all_results <- bind_rows(rol_choice, maxdiff_choice, rol_choice_no_ties)

ggplot(all_results, aes(x = True_value_normalized, y = median, color = Dataset)) +
  geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_point(aes(y = median), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Normalized True Utility",
       y = "Normalized Estimated Utility",
       title = "Normalized Utility Estimates vs True Values",
       subtitle = "With interior 90% credibility intervals") +
  scale_color_manual(values = c("Rank-Ordered Logit w ties" = "green", 
                                "MaxDiff" = "blue", 
                                "Rank-Ordered Logit (no ties)" = "red")) +
  theme_minimal() +
  coord_fixed(ratio = 1) +
  facet_wrap(~Dataset)

rmse_rol <- sqrt(mean((rol_choice$True_value_normalized - rol_choice$median)^2))
rmse_maxdiff <- sqrt(mean((maxdiff_choice$True_value_normalized - maxdiff_choice$median)^2))
rmse_rol_no_ties <- sqrt(mean((rol_choice_no_ties$True_value_normalized - rol_choice_no_ties$median)^2))

print(paste("RMSE Rank-Ordered Logit with ties:", round(rmse_rol, 4)))
print(paste("RMSE MaxDiff:", round(rmse_maxdiff, 4)))
print(paste("RMSE Rank-Ordered Logit (no ties):", round(rmse_rol_no_ties, 4)))

plot_with_rmse <- last_plot() +
  labs(subtitle = paste("With interior 90% credibility intervals\n",
                        "RMSE ROL w/ ties:", round(rmse_rol, 4), "\n",
                        "RMSE MaxDiff:", round(rmse_maxdiff, 4), "\n",
                        "RMSE ROL (no ties):", round(rmse_rol_no_ties, 4)))

print(plot_with_rmse)
