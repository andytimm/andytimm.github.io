library(tidyverse)
library(loo)

options(mc.cores = 8)

simulated_data <- tibble(
q_x = rnorm(100000),
manageable_p_x = rt(100000,5),
unmanageable_p_x = rcauchy(100000),
manageable_ratios = dt(q_x,5)/dnorm(q_x),
unmanageable_ratios = dcauchy(q_x)/dnorm(q_x)
)

mean(simulated_data$manageable_ratios)
max(simulated_data$manageable_ratios)

mean(simulated_data$unmanageable_ratios)
max(simulated_data$unmanageable_ratios)

reasonable_ratios <- simulated_data %>%
  arrange(manageable_ratios) %>%
  mutate(n = seq(1,100000)) %>%
  ggplot(aes(x = n,y = manageable_ratios)) +
  geom_point() +
  ggtitle("A pretty typical 'close to good' set of importance ratios")

unreasonable_ratios <- simulated_data %>%
  arrange(unmanageable_ratios) %>%
  mutate(n = seq(1,100000)) %>%
  ggplot(aes(x = n,y = unmanageable_ratios)) +
  geom_point() +
  ggtitle("A pretty typical 'unsaveable' set of importance ratios")

density_visual <- simulated_data %>%
                    pivot_longer(c(q_x,manageable_p_x,unmanageable_p_x),
                                 values_to = "draws",
                                 names_to = "distributions") %>%
                    ggplot(aes(x = draws, color = distributions)) +
                    geom_density() +
                    xlim(-10,10) +
                    ggtitle("Visualizing the distributions in question") +
                    theme(legend.position="none")

manageable_psis <- psis(log(simulated_data$manageable_ratios),
                       r_eff = NA)

manageable_psis$diagnostics

unmanageable_psis <- psis(log(simulated_data$unmanageable_ratios),
                          r_eff = NA)

unmanageable_psis$diagnostics

