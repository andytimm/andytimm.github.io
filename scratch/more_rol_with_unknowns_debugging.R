library(tidyverse)
library(rstan)

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
ranked_options <- crossing(individual = 1:I, task = 1:Tasks, option = 1:J) %>% 
  mutate(row = 1:n()) %>%
  group_by(individual, task) %>% 
  mutate(
    fixed_utility = as.numeric(X[row,] %*% as.numeric(beta_i[first(individual),])),
    plus_gumbel_error = fixed_utility + rgumbel(n()),
    true_rank = rank(-plus_gumbel_error),
    true_order = order(-plus_gumbel_error),
    # Simulate partial ordering by only observing top 1 and bottom 1
    observed_order = case_when(
      true_order <= 1 ~ true_order,
      true_order == J ~ J,
      TRUE ~ 3  # tie all middle orders
    ),
    # For original models
    best_choice = as.numeric(true_order == 1),
    worst_choice = as.numeric(true_order == J)
  ) %>%
  ungroup()

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
    vector[K] exp_delta = exp(delta - max(delta));  // Subtract max for numerical stability
    real out = 0;
    int current_pos = 1;
    real remaining_sum = sum(exp_delta);
    real epsilon = 1e-10;  // Small constant to avoid division by zero

    //print(\"Initial exp_delta: \", exp_delta);
    //print(\"Initial remaining_sum: \", remaining_sum);

    while (current_pos <= K) {
      int tied_count = 1;
      real tied_sum = exp_delta[rank_order[current_pos]];
      
      while (current_pos + tied_count <= K && rank_order[current_pos] == rank_order[current_pos + tied_count]) {
        tied_sum += exp_delta[rank_order[current_pos + tied_count]];
        tied_count += 1;
      }

      //print(\"Current position: \", current_pos);
      //print(\"Tied count: \", tied_count);
      //print(\"Tied sum: \", tied_sum);

      if (tied_count == 1) {
        out += delta[rank_order[current_pos]] - log(fmax(remaining_sum, epsilon));
      } else {
        out += log(fmax(tied_sum, epsilon)) - log(fmax(remaining_sum, epsilon));
      }

      //print(\"Log-likelihood contribution: \", out);

      remaining_sum = fmax(remaining_sum - tied_sum, 0.0);  // Ensure non-negativity
      current_pos += tied_count;

      //print(\"Updated remaining_sum: \", remaining_sum);
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
  // In R, this is observed_order within each task
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
  matrix[I, P] beta_individual = rep_matrix(beta\', I) + X2 * Gamma\' + z * diag_pre_multiply(tau, L_Omega);
}
model {
  tau ~ normal(0, .5);
  beta ~ normal(0, 1);
  to_vector(z) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(Gamma) ~ normal(0, 1);

  for(t in 1:T) {
    vector[K] utilities;
    utilities = X[start[t]:end[t]] * beta_individual[task_individual[t]]\';

    // Center utilities for numerical stability
    utilities = utilities - mean(utilities);

    // Scale utilities
    utilities = utilities / max(fabs(utilities));

    // Print utilities and their transformations for debugging
    //print(\"Task: \", t);
    //print(\"Original utilities: \", X[start[t]:end[t]] * beta_individual[task_individual[t]]\');
    //print(\"Centered utilities: \", utilities + mean(utilities));
    //print(\"Scaled utilities: \", utilities);
    //print(\"Rank order: \", rank_order[start[t]:end[t]]);

    // Use the improved likelihood function
    real ll = rank_logit_ties_lpmf(rank_order[start[t]:end[t]] | utilities);
    //print(\"Log-likelihood: \", ll);
    if (!is_nan(ll) && !is_inf(ll)) {
      target += ll;
    } else {
      //print(\"Warning: Invalid log-likelihood at task \", t);
    }
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

# Print a summary of the results
print(fit, pars = c("beta", "tau"))

# Extract the posterior samples
posterior_samples <- extract(fit)

# Calculate mean betas and 95% credible intervals
beta_summary <- data.frame(
  mean_beta = apply(posterior_samples$beta, 2, mean),
  lower_ci = apply(posterior_samples$beta, 2, quantile, probs = 0.025),
  upper_ci = apply(posterior_samples$beta, 2, quantile, probs = 0.975)
)

print(beta_summary)

# If you want to examine individual-level betas:
beta_individual_summary <- apply(posterior_samples$beta_individual, c(2, 3), mean)
print(head(beta_individual_summary))

# Assuming you have the true beta_i values stored in a matrix called 'beta_i'
# and the estimated values are in the 'fit' object from Stan

# Extract the posterior means for beta_individual
posterior_means <- summary(fit)$summary[grep("beta_individual", rownames(summary(fit)$summary)), "mean"]

# Reshape the posterior means into a matrix
estimated_beta_i <- matrix(posterior_means, nrow = I, ncol = P, byrow = TRUE)

# Calculate the difference between estimated and true values
differences <- estimated_beta_i - beta_i

# Calculate RMSE for each parameter
rmse <- sqrt(colMeans(differences^2))

# Calculate correlation between estimated and true values for each parameter
correlations <- sapply(1:P, function(j) cor(estimated_beta_i[,j], beta_i[,j]))

# Print results
cat("RMSE for each parameter:\n")
print(rmse)

cat("\nCorrelation between estimated and true values for each parameter:\n")
print(correlations)

# Optionally, create a plot to visualize the comparison
library(ggplot2)

plot_data <- data.frame(
  True = as.vector(beta_i),
  Estimated = as.vector(estimated_beta_i),
  Parameter = rep(paste("Parameter", 1:P), each = I)
)

ggplot(plot_data, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  facet_wrap(~ Parameter, scales = "free") +
  labs(title = "Comparison of True vs Estimated Utilities",
       x = "True Utility",
       y = "Estimated Utility")
