library(tidyverse)

rules_1 <- rules_1[order(rules_1$time), ]

rules_1

n_pop <- 1e4

time <- 1:10

niterations <- 1000

set.seed(4730)

### Creating objects to store results
mylist <- list()
listprev <- list()

mean_e <- rep(NA, niterations)
mean_d <- rep(NA, niterations)
mean_c <- rep(NA, niterations)
RR_EC <- rep(NA, niterations)
RR_DC <- rep(NA, niterations)
RR_DEC <- rep(NA, niterations)

for(i in 1:niterations) {
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
  for(j in unique(df$time)){
    df$u1[df$time == j] <- rbinom(n_pop, 1, rules_1$u1[j])
    df$u2[df$time == j] <- rbinom(n_pop, 1, rules_1$u2[j])
    df$u3[df$time == j] <- rbinom(n_pop, 1, rules_1$u3[j])
    df$u4[df$time == j] <- rbinom(n_pop, 1, rules_1$u4[j])
  }
  
  ### Simulating E = 1 if U4 OR C and U3, and D = 1 if C and U1 OR U2
  df <- df %>%
    mutate(e = ifelse((u4 == 1) | (c == 1 & u3 == 1), 1, 0), 
           d = ifelse((c == 1 & u1 == 1) | (u2 == 1), 1, 0))
  
  mean_e[i] <- mean(df$e)
  mean_d[i] <- mean(df$d)
  mean_c[i] <- mean(df$c)
  
  listprev[[i]] <- df %>% group_by(time) %>%
                          summarise(pr_c = mean(c), 
                                    pr_e = mean(e), 
                                    pr_d = mean(d))
}

mean_e

### Function to select C, E, or D from each df in the list of dataframes
list_select_pr <- function(some_list, x, y){
  return(some_list[[x]][[y]])
}

lapply(mylist, aggregate, list(df$time), mean)

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

