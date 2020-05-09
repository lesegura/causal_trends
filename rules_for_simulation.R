##################################################################
#                                                                #
#                     Rules for Simulation                       #
#                         Scenario 1                             #
#                                                                #
# Rule 1: The prevalence of C is constant at each time t         #
# Rule 2: C causes D and E                                       #
# Rule 3: The prevalence of U1, & U2 are constant over time      #
# Rule 4: We are assuming monotonicity for causes of E & D       #
# Rule 5: P4 is constant for E at 0.60 and for D at 0.70         #
# Rule 6: The prevalence of the totality of the causes of E      #
#         is 1 - P4 for E. The prevalence of all causes of D     #
#         is 1 - P4 for D                                        #
# Rule 6: The prevalence of U4 increases over time               #
# Rule 7: The prevalence of U3 is stable over time               #
# Rule 8: Pr(E) = 1 if U3 or C and U4                            #
# Rule 9: Pr(D) = 1 if C and U1 or U2. E does not cause D        #
#                                                                #
##################################################################


library(tidyverse)

### Exponential growth function
exp_fun <- function(x, r, t){
  x * ((1 + r) ^ t)
}


### Prevalence of Causal Partners

u4 <- rep(NA, 10)

for(i in 0:9){
  u4[i + 1] <- exp_fun(0.04, 0.285, i)
}

rules_1 <- tibble(time = 1:10, 
                u1 = 0.25, 
                u2 = 0.05, 
                ### substract the middle of the venn diagram (double counted)
                RR_D_C = (u1 + u2 - (u1*u2)) / (u2 - (u1*u2)),
                u3 = 0.4 - u4, 
                u4 = u4, 
                RR_E_C = (u3 + u4 - (u3*u4)) / (u4 - (u3*u4)))

rules_1

##################################################################
#                                                                #
#                     Rules for Simulation                       #
#                         Scenario 2                             #
#                                                                #
# Rule 1: The prevalence of C is constant at each time t         #
# Rule 2: C causes D and E                                       #
# Rule 3: The prevalence of U3, & U4 are constant over time      #
# Rule 4: We are assuming monotonicity for causes of E & D       #
# Rule 5: P4 is constant for E at 0.60 and for D at 0.70         #
# Rule 6: The prevalence of the totality of the causes of E      #
#         is 1 - P4 for E. The prevalence of all causes of D     #
#         is 1 - P4 for D                                        #
# Rule 6: The prevalence of U1 decreases over time               #
# Rule 7: The prevalence of U2 increases over time               #
# Rule 8: Pr(E) = 1 if U3 or C and U4                            #
# Rule 9: Pr(D) = 1 if C and U1 or U2. E does not cause D        #
#                                                                #
##################################################################

### Prevalence of Causal Partners

u2 <- rep(NA, 10)

for(i in 0:9){
  if(i > 7 ) {
    u2[i + 1] <- exp_fun(0.05, 0.285, 7)
  } else {
    u2[i + 1] <- exp_fun(0.05, 0.285, i)
  }
}

rules_2 <- tibble(time = 1:10, 
                  u1 = 0.30 - u2, 
                  u2 = u2, 
                  RR_D_C = (u1 + u2) / u2,
                  u3 = 0.36, 
                  u4 = 0.40 - u3, 
                  RR_E_C = (u3 + u4) / u4)

rules_2

##################################################################
#                                                                #
#                     Rules for Simulation                       #
#                         Scenario 3                             #
#                                                                #
# Rule 1: The prevalence of C is constant at each time t         #
# Rule 2: C causes D and E                                       #
# Rule 3: We are assuming monotonicity for causes of E & D       #
# Rule 4: P4 is constant for E at 0.60 and for D at 0.70         #
# Rule 5: The prevalence of the totality of the causes of E      #
#         is 1 - P4 for E. The prevalence of all causes of D     #
#         is 1 - P4 for D                                        #
# Rule 6: The prevalence of U1 & U3 decreases over time          #
# Rule 7: The prevalence of U2 & U4 increases over time          #
# Rule 8: Pr(E) = 1 if U3 or C and U4                            #
# Rule 9: Pr(D) = 1 if C and U1 or U2. E does not cause D        #
#                                                                #
##################################################################

### Prevalence of Causal Partners

### Prevalence of Causal Partners

u2 <- rep(NA, 10)

for(i in 0:9){
  if(i > 7 ) {
    u2[i + 1] <- exp_fun(0.05, 0.285, 7)
  } else {
    u2[i + 1] <- exp_fun(0.05, 0.285, i)
  }
}

u4 <- rep(NA, 10)

for(i in 0:9){
  u4[i + 1] <- exp_fun(0.04, 0.285, i)
}

rules_3 <- tibble(time = 1:10, 
                  u1 = 0.30 - u2, 
                  u2 = u2, 
                  RR_D_C = (u1 + u2) / u2,
                  u3 = 0.4 - u4, 
                  u4 = u4, 
                  RR_E_C = (u3 + u4) / u4)

rules_3
