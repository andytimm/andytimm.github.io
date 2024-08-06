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
  J = J,
  N = I * Tasks,  # total number of ranking tasks
  K = P,  # number of covariates
  orders = ranked_options %>%
    select(individual, task, option, observed_order) %>%
    pivot_wider(names_from = option, values_from = observed_order) %>%
    select(-individual, -task) %>%
    as.matrix(),
  X = X
)


efron_simplified_rol <- "
data {
  int<lower=2> J;  // number of items to order
  int<lower=1> N;  // number of orderings
  int<lower=1> K;  // number of covariates
  array[N, J] int<lower=1, upper=J> orders;  // Partial orders
  matrix[N*J, K] X;  // covariate matrix
}

parameters {
  vector[K] beta;  // regression coefficients
}

model {
  beta ~ normal(0, 1);  // prior

  for (n in 1:N) {
    vector[J] mu;
    for (j in 1:J) {
      mu[j] = X[(n-1)*J + j] * beta;
    }
    
    // Center mu to improve numerical stability
    mu = mu - mean(mu);
    
    int current_order = 1;
    real log_lik = 0;
    
    while (current_order <= J) {
      int tied_count = 1;
      vector[J] exp_mu = exp(mu);
      real tied_sum = exp_mu[orders[n, current_order]];
      
      // Find ties
      while (current_order + tied_count <= J && 
             orders[n, current_order] == orders[n, current_order + tied_count]) {
        tied_sum += exp_mu[orders[n, current_order + tied_count]];
        tied_count += 1;
      }
      
      real remaining_sum = sum(exp_mu) - sum(exp_mu[orders[n, 1:(current_order-1)]]);
      
      if (tied_count == 1) {
        // No tie, standard logit term
        log_lik += mu[orders[n, current_order]] - log(remaining_sum);
      } else {
        // Tied (unknown) orderings
        log_lik += log(tied_sum) - log(remaining_sum);
      }
      
      current_order += tied_count;
    }
    
    target += log_lik;
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

# Print a summary of the results
print(fit, pars = "utility")

# Extract the posterior samples
posterior_samples <- extract(fit)

# Calculate mean utilities and 95% credible intervals
utility_summary <- data.frame(
  mean_utility = apply(posterior_samples$utility, 2, mean),
  lower_ci = apply(posterior_samples$utility, 2, quantile, probs = 0.025),
  upper_ci = apply(posterior_samples$utility, 2, quantile, probs = 0.975)
)

print(utility_summary)

# Debugging

# 1. Check convergence diagnostics
print(fit, pars = c("utility"))

# Plot trace plots - seems fine
rstan::traceplot(fit, pars = c("utility"))

# 1. Extract true utilities
true_utilities <- colMeans(beta_i)

# 2. Calculate estimated utilities
estimated_utilities <- summary(fit, pars = "utility")$summary[, "mean"]

# 3. Compare true vs. estimated utilities
plot(true_utilities, estimated_utilities, 
     xlab = "True Utilities", ylab = "Estimated Utilities", 
     main = "Estimated vs True Utilities",
     pch = 19, col = "blue") +
abline(a = 0, b = 1, col = "red", lty = 2) +  # Add y=x line
text(true_utilities, estimated_utilities, labels = 1:length(true_utilities), pos = 4)

print(str(stan_data_list))
print(summary(stan_data_list))

fit <- sampling(compiled_model, 
                data = stan_data_list, 
                chains = 1, 
                iter = 10,  # Small number of iterations for quick debugging
                warmup = 5,
                algorithm = "Fixed_param")  # This algorithm just draws from the prior


# Now try full Efron

efron_full_rol <- "
data {
  int<lower=2> K;  // Number of choice options
  int<lower=1> N;  // Number of choice tasks
  int<lower=1> I;  // Number of individuals
  int<lower=1> P;  // Number of covariates
  array[N] int<lower=1, upper=K> best;  // Index of best choice for each task
  array[N] int<lower=1, upper=K> worst;  // Index of worst choice for each task
  array[N] int<lower=1, upper=I> individual;  // Individual for each task
  matrix[N, P] X;  // Covariate matrix for each choice task
}

parameters {
  matrix[K-1, P] beta;  // Coefficients for covariates
  vector[K-1] alpha;  // Intercepts for utilities
  matrix[I, K-1] z;  // Individual-level random effects
  vector<lower=0>[K-1] tau;  // Standard deviations of random effects
  cholesky_factor_corr[K-1] L_Omega;  // Cholesky factor of correlation matrix
}

transformed parameters {
  matrix[I, K] individual_utility;
  
  for (i in 1:I) {
    individual_utility[i, 1] = 0;  // Reference level
    individual_utility[i, 2:K] = alpha + (tau .* (L_Omega * z[i]'))';
  }
}

model {
  // Priors
  to_vector(beta) ~ normal(0, 1);
  alpha ~ normal(2, 5);
  to_vector(z) ~ normal(0, 1);
  tau ~ cauchy(0, 2.5);
  L_Omega ~ lkj_corr_cholesky(2);
  
  // Likelihood
  for (n in 1:N) {
    for (n in 1:N) {
    int i = individual[n];
    vector[K] task_utility;
    
    for (k in 1:K) {
      task_utility[k] = individual_utility[i, k] + dot_product(X[n], beta[k,]);
    }
    
    real log_lik = 0;
    
    // Contribution from best choice
    log_lik += task_utility[best[n]] - log_sum_exp(task_utility);
    
    // Contribution from worst choice
    vector[K-1] remaining_utility;
    if (worst[n] == 1) {
      remaining_utility = task_utility[2:K];
    } else if (worst[n] == K) {
      remaining_utility = task_utility[1:(K-1)];
    } else {
      remaining_utility = append_row(
        task_utility[1:(worst[n]-1)], 
        task_utility[(worst[n]+1):K]
      );
    }
    log_lik += -task_utility[worst[n]] - log_sum_exp(-remaining_utility);
    
    // Efron's approximation for middle choices
    int n_middle = K - 2;
    if (n_middle > 0) {
      vector[n_middle] middle_utilities;
      int pos = 1;
      for (k in 1:K) {
        if (k != best[n] && k != worst[n]) {
          middle_utilities[pos] = task_utility[k];
          pos += 1;
        }
      }
      
      real sum_middle = sum(middle_utilities);
      for (j in 1:n_middle) {
        log_lik += (sum_middle - j * log(j)) / n_middle;
      }
    }
    
    target += log_lik;
  }
}

generated quantities {
  matrix[K, P] beta_full;
  matrix[I, K] utility_centered;
  
  beta_full[1] = rep_row_vector(0, P);
  beta_full[2:K] = beta;
  
  for (i in 1:I) {
    utility_centered[i] = individual_utility[i] - mean(individual_utility[i]);
  }
}
"

ranked_options <- indexes %>% 
  group_by(individual, task) %>% 
  mutate(
    fixed_utility = as.numeric(X[row,] %*% as.numeric(beta_i[first(individual),])),
    plus_gumbel_error = fixed_utility + rgumbel(n()),
    rank = rank(-plus_gumbel_error),
    order = order(rank),
    best = which.max(plus_gumbel_error),
    worst = which.min(plus_gumbel_error)
  ) %>%
  ungroup()

# Prepare data for Stan
stan_data <- list(
  K = J,  # Number of choice options
  N = nrow(ranked_options),  # Total number of choice tasks
  I = max(ranked_options$individual),  # Number of individuals
  P = ncol(X),  # Number of covariates
  best = ranked_options$best,
  worst = ranked_options$worst,
  individual = ranked_options$individual,
  X = X  # Your covariate matrix
)

# Compile and fit the model
fit <- sampling(stan_model(model_code = efron_full_rol), data = stan_data, 
                chains = 4, iter = 2000, warmup = 1000)

# Extract and summarize results
print(fit, pars = c("alpha", "beta", "tau"))