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