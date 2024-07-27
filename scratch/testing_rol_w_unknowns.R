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

# Function to generate all permutations of middle options
generate_permutations <- function(n) {
  if (n == 2) return(matrix(c(1,2, 2,1), ncol=2, byrow=TRUE))
  perms <- combinat::permn(1:n)
  do.call(rbind, perms)
}

# Generate choice sets with best, worst, and unknown middle options
choice_sets <- ranked_options %>%
  group_by(individual, task) %>%
  mutate(
    is_best = rank == 1,
    is_worst = rank == max(rank),
    is_middle = !is_best & !is_worst
  ) %>%
  summarise(
    best = which(is_best),
    worst = which(is_worst),
    middle = list(which(is_middle)),
    .groups = 'drop'
  )

# Generate all permutations of middle options
middle_perms <- generate_permutations(J - 2)

# Prepare data for Stan
stan_data <- list(
  N = nrow(X),
  T = nrow(choice_sets),
  I = I,
  P = P,
  P2 = P2,
  J = J,
  K = J - 2,  # number of unknown middle options
  X = X,
  X2 = W,
  best = choice_sets$best,
  worst = choice_sets$worst,
  middle_perms = middle_perms,
  task = tt$task_number,
  task_individual = tt$individual,
  start = tt$start,
  end = tt$end
)

rol_w_unknowns <- "
functions {
  // tgamma(n) = factorial(n-1), so this should work fine
  int factorial(int n) {
    return tgamma(n + 1);
  }

  real partial_rank_log_prob(array[] int ranks, vector utilities, array[,] int perms) {
    int n_perms = size(perms);
    vector[n_perms] perm_probs;
    int J = size(ranks);
    
    for (p in 1:n_perms) {
      vector[J] perm_utilities;
      real prob = 0;
      perm_utilities[1] = utilities[ranks[1]];  // best
      perm_utilities[J] = utilities[ranks[J]];  // worst
      for (i in 2:(J-1)) {
        perm_utilities[i] = utilities[ranks[perms[p, i-1]]];
      }
      
      for (j in 1:(J-1)) {
        prob += perm_utilities[j] - log_sum_exp(perm_utilities[j:J]);
      }
      perm_probs[p] = prob;
    }
    return log_sum_exp(perm_probs);
  }
}

data {
  int<lower=2> N;  // number of observations
  int<lower=1> T;  // number of tasks
  int<lower=1> I;  // number of individuals
  int<lower=1> P;  // number of predictors
  int<lower=1> P2; // number of individual-level predictors
  int<lower=3> J;  // number of alternatives per task
  int<lower=1> K;  // number of unknown middle options
  matrix[N, P] X;  // predictors
  matrix[I, P2] X2; // individual-level predictors
  array[T] int<lower=1,upper=J> best;  // index of best option for each task
  array[T] int<lower=1,upper=J> worst;  // index of worst option for each task
  array[factorial(K), K] int<lower=1,upper=J-2> middle_perms;  // all permutations of middle options
  array[T] int<lower=1,upper=T> task;  // task index
  array[T] int<lower=1,upper=I> task_individual;  // individual index for each task
  array[T] int<lower=1,upper=N> start;  // start index for each task
  array[T] int<lower=1,upper=N> end;  // end index for each task
}

parameters {
  vector[P] beta;
  matrix[P, P2] Gamma;
  vector<lower=0>[P] tau;
  matrix[I, P] z;
  cholesky_factor_corr[P] L_Omega;
}

transformed parameters {
  matrix[I, P] beta_individual = rep_matrix(beta', I) + X2 * Gamma' + z * diag_pre_multiply(tau, L_Omega);
}

model {
  tau ~ normal(0, .5);
  beta ~ normal(0, .5);
  to_vector(z) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(Gamma) ~ normal(0, 1);

  for (t in 1:T) {
    vector[J] utilities;
    array[J] int ranks;
    utilities = X[start[t]:end[t]] * beta_individual[task_individual[t]]';
    ranks[1] = best[t];
    ranks[J] = worst[t];
    int idx = 2;
    for (j in 1:J) {
      if (j != best[t] && j != worst[t]) {
        ranks[idx] = j;
        idx += 1;
      }
    }
    target += partial_rank_log_prob(ranks, utilities, middle_perms);
  }
}
"

compiled_rol_unk_model <- stan_model(model_code = rol_w_unknowns)
rol_fit <- sampling(compiled_rol_unk_model,
                    data = stan_data,
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
