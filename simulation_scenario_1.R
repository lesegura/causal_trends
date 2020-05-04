library(tidyverse)

rules_1 <- rules_1[order(rules_1$time), ]

rules_1

n_pop <- 1e4

time <- 1:10

niterations <- 1000

set.seed(4730)

### Creating storing objects
mean_e <- rep()

### Creating the dataset with 10,000 observations per time. 
df <- tibble(time = rep(time, each = n_pop), 
             ### Prevalence of C is time invariant at 20%
               c = rbinom(n_pop * time, 1, .2), 
               u1 = NA, 
               u2 = NA, 
               u3 = NA, 
               u4 = NA)

### Sampling the causal partners from a binomial distribution with time varying probabilities
### according to rules_1
for(i in 1:10){
    df$u1[df$time == i] <- rbinom(n_pop, 1, rules_1$u1[i])
    df$u2[df$time == i] <- rbinom(n_pop, 1, rules_1$u2[i])
    df$u3[df$time == i] <- rbinom(n_pop, 1, rules_1$u3[i])
    df$u4[df$time == i] <- rbinom(n_pop, 1, rules_1$u4[i])
}

### Simulating E = 1 if U4 OR C and U3, and D = 1 if C and U1 OR U2
df <- df %>%
    mutate(e = ifelse((u4 == 1) | (c == 1 & u3 == 1), 1, 0), 
           d = ifelse((c == 1 & u1 == 1) | (u2 == 1), 1, 0))


df %>%
  summarise(pr_c = mean(c),
            pr_e = mean(e),
            pr_d = mean(d))

df %>%
  group_by(time) %>%
  summarise(pr_u1 = mean(u1),
            pr_u2 = mean(u2), 
            pr_u3 = mean(u3), 
            pr_u4 = mean(u4), 
            pr_c = mean(c),
            pr_e = mean(e),
            pr_d = mean(d))

  




rm(mylist)

sim_data_1$u3 <- NA
sim_data_1$u4 <- NA

for(i in 1:10){
  sim_data_1$u3[sim_data_1$time == i] <- rbinom(n_pop * time, 1, rules_1$u3[i])
}

for(i in 1:10){
  sim_data_1$u4[sim_data_1$time == i] <- rbinom(n_pop * time, 1, rules_1$u4[i])
}

sim_data_1 <- sim_data_1 %>%
  mutate(e = ifelse((u3 == 1) | (c == 1 & u4 == 1), 1, 0), 
         d = ifelse((c == 1 & u1 == 1) | (u2 == 1), 1, 0))

sim_data_1 %>%
  group_by(time) %>%
  summarise(pr_u1 = mean(u1),
            pr_u2 = mean(u2), 
            pr_u3 = mean(u3), 
            pr_u4 = mean(u4), 
            pr_c = mean(c),
            pr_e = mean(e),
            pr_d = mean(d))

rules_1

