rm(list=ls())

library(formatR)
library(knitr)
library(ggplot2)
library(dagitty)
library(pander)

setwd("/Users/luissegura/Dropbox/Causal Trends")

### FUNCTIONS
# RR
risk_r<-function(p1, p2){
  (p1+p2) / p1
}

# Odds Function
odds<-function(x){
  x/(1-x)
}

# OR
or_<-function(p1, p2){
  ((p1+p2) / (p1)) / ((1-(p1+p2)) / ((1-(p1+p2)) + p2))
}

# Log Odds Function
ln_odds<-function(x){
  log(x/(1-x))
}

# Exponential growth function
exp_fun<-function(x, r, t){
  x*((1+r)^t)
}

# Exponentiated Logit Function (calculates a probability from an odds ratio)
logit<-function(a, b, x){
  exp(a + b*x) / (1 + exp(a + b*x))
}


######################### RULES FOR SIMULATION ############################

#---------------------------------------------#
#                                             #
#                 CAUSES OF E                 #
#                                             #
#---------------------------------------------#
cause_tab2<-data.frame(matrix(NA, nrow=10, ncol=1))
cause_tab2[, 1]<-2009:2018

### Arbitrary selection of the pr of p1s (0.04) and p2s (0.36)
cause_tab2[, 2]<- rep(c(.04))
cause_tab2[, 3]<- rep(c(.36))
cause_tab2[, 4]<- risk_r(.04, 0.36)
cause_tab2[, 5]<- or_(.04, 0.36)

#---------------------------------------------#
#                                             #
#                 CAUSES OF D                 #
#                                             #
#---------------------------------------------#

p1_d<-data.frame(matrix(NA, nrow=10, ncol=1))
p1_d[1,]<-exp_fun(0.05, 0.285, 0)
p1_d[2,]<-exp_fun(0.05, 0.285, 1)
p1_d[3,]<-exp_fun(0.05, 0.285, 2)
p1_d[4,]<-exp_fun(0.05, 0.285, 3)
p1_d[5,]<-exp_fun(0.05, 0.285, 4)
p1_d[6,]<-exp_fun(0.05, 0.285, 5)
p1_d[7,]<-exp_fun(0.05, 0.285, 6)
p1_d[8,]<-exp_fun(0.05, 0.285, 7)
p1_d[9,]<-0.2892623
p1_d[10,]<-0.2892623


cause_tab2[, 6]<- p1_d

### Sequence of Prevalences of p2s
### Note here 0.7 hard coded as the prevalence of p4s
p2_d<-.3-p1_d
cause_tab2[, 7]<- p2_d

### Estimating the RR and OR
cause_tab2[, 8]<- rr_d<-risk_r(p1_d, p2_d)
cause_tab2[, 9]<- or_d<-or_(p1_d, p2_d)

colnames(cause_tab2)<-c("year", "p1_for_E", "p2_for_E", "RR_E", "OR_E", "p1_for_D", "p2_for_D", "RR_D", "OR_D")

pander(cause_tab2, caption = "Causes (p1s and p2s) of E and D at each T=t and their expected magnitude of effect")

#---------------------------------------------#
#                                             #
#           SETTING THE PARAMETERS            #
#                                             #
#---------------------------------------------#


betas_gammas2<- data.frame(cause_tab2$year)
colnames(betas_gammas2)<-"year"
betas_gammas2$g0_for_E<- ln_odds(cause_tab2$p1_for_E)
betas_gammas2$g1_for_E<- log(cause_tab2$OR_E)
betas_gammas2$b0_for_D<- ln_odds(cause_tab2$p1_for_D)
betas_gammas2$b1_for_D<- log(cause_tab2$OR_D)

pander(betas_gammas2, caption = "Parameters that we are going to use for the second simulation")


#---------------------------------------------------------#
#                                                         #
#          SIMULATING DATASET (100 ITERATIONS)            #
#                                                         #
#---------------------------------------------------------#

### Creating the dataframe of all 100 datasets
listofdf<-list() ### list of dataframes

### Looping the creation of 100 datasets 
for(x in 1:100){
  df<-data.frame(rep(c(2009:2018), each=1000)) ## this variable is an indicator of calendar year, 2008 - 2018.
  colnames(df)<-"year"
  ### Simulating C
  df$confounder<-rbinom(nrow(df), 1, .2)
  ### Estimating the Probability of Exposure
  for(i in unique(df$year)){
    df$p_exposure[df$year==i]<-logit(betas_gammas2$g0_for_E[betas_gammas2$year==i], 
                                     betas_gammas2$g1_for_E[betas_gammas2$year==i], 
                                     df$confounder[df$year==i])
  }
  ### Simulating E
  df$exposure<-rbinom(nrow(df), 1, pr=df$p_exposure)
  ### Estimating Probability of D
  for(i in unique(df$year)){
    df$p_disease[df$year==i]<-logit(betas_gammas2$b0_for_D[betas_gammas2$year==i], 
                                    betas_gammas2$b1_for_D[betas_gammas2$year==i], 
                                    df$confounder[df$year==i])
  }
  ## Simulating D
  df$disease<-rbinom(nrow(df), 1, pr=df$p_disease)
  
  listofdf[[x]]<-df
}


### Naming the dataframes in the list
for(i in 1:length(listofdf)){
  names(listofdf)[i]<-paste("df", i, sep = "")
}


### Function to select C, E, or D from each df in the list of dataframes
list_select_pr<-function(some_list, x, y){
  return(some_list[[x]][[y]])
}

for(i in 1:100){
  c_all<-mean(list_select_pr(listofdf, i, "confounder"))
  e_all<-mean(list_select_pr(listofdf, i, "exposure"))
  d_all<-mean(list_select_pr(listofdf, i, "disease"))
}

tab_pr<-data.frame(cbind(c_all, e_all, d_all))
pander(tab_pr, caption = "Prevalence of C, E, and D after 100 iterations of the simulated dataset")


#------------------------------------------------------#
#                                                      #
#         Prevalence of C, E, and D by year            #
#                                                      #
#------------------------------------------------------#


### Generating a table of the mean of C, E, and D for each simulated dataset by year
listofpr<-lapply(listofdf, aggregate, list(df$year), mean)

### A dataframe of the confounder in all 100 simulations
df.c<-data.frame(matrix(NA, nrow=10, ncol=100))

for(i in 1:100){
  df.c[, i]<-listofpr[[i]][["confounder"]]
  colnames(df.c)[i]<-paste("sim", i, sep = "")
}

df.c[, 101]<-2009:2018
colnames(df.c)[101]<-"year"
df.c<-df.c[, c(101, 1:100)]
df.c<-data.frame(t(df.c[, -1]))
colnames(df.c)<-c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10")


### Grand mean of C from 100 simulated datasets at each year 
df.master<-data.frame(c(2009:2018))
df.master[, 2]<-data.frame(rbind(mean(df.c$t1), mean(df.c$t2), mean(df.c$t3), mean(df.c$t4), mean(df.c$t5), 
                                 mean(df.c$t6), mean(df.c$t7), mean(df.c$t8), mean(df.c$t9), mean(df.c$t10)))

colnames(df.master)<-c("year", "confounder")

### A dataframe of E in all 100 simulations
df.e<-data.frame(matrix(NA, nrow=10, ncol=100))

for(i in 1:100){
  df.e[, i]<-listofpr[[i]][["exposure"]]
  colnames(df.e)[i]<-paste("sim", i, sep = "")
}

df.e[, 101]<-2009:2018
colnames(df.e)[101]<-"year"
df.e<-df.e[, c(101, 1:100)]
df.e<-data.frame(t(df.e[, -c(1)]))
colnames(df.e)<-c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10")

### Grand mean of C from 100 simulated datasets at each year 
df.master[, 3]<-data.frame(rbind(mean(df.e$t1), mean(df.e$t2), mean(df.e$t3), mean(df.e$t4), mean(df.e$t5), 
                                 mean(df.e$t6), mean(df.e$t7), mean(df.e$t8), mean(df.e$t9), mean(df.e$t10)))

colnames(df.master)[3]<-c("exposure")
scenario2.pr<-df.master


### A dataframe of D in all 100 simulations
df.d<-data.frame(matrix(NA, nrow=10, ncol=100))

for(i in 1:100){
  df.d[, i]<-listofpr[[i]][["disease"]]
  colnames(df.d)[i]<-paste("sim", i, sep = "")
}

df.d[, 101]<-2009:2018
colnames(df.d)[101]<-"year"
df.d<-df.d[, c(101, 1:100)]
df.d<-data.frame(t(df.d[, -c(1)]))
colnames(df.d)<-c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10")

### Grand mean of D from 100 simulated datasets at each year 
df.master[, 4]<-data.frame(rbind(mean(df.d$t1), mean(df.d$t2), mean(df.d$t3), mean(df.d$t4), mean(df.d$t5), 
                                 mean(df.d$t6), mean(df.d$t7), mean(df.d$t8), mean(df.d$t9), mean(df.d$t10)))

colnames(df.master)[4]<-c("disease")

pander(df.master, caption = "Prevalence of C, E, and D from 100 simulated datasets")

df.plot<-data.frame(matrix(NA, nrow=30, ncol=3))
df.plot[, 1]<-df.master$year
df.plot[, 2]<-rep(c("confounder", "exposure", "disease"), each = 10)
df.plot[1:10, 3]<-df.master$confounder
df.plot[11:20, 3]<-df.master$exposure
df.plot[21:30, 3]<-df.master$disease
colnames(df.plot)<-c("year", "variable", "prob")

ggplot(data = df.plot, 
       aes(x = as.numeric(year),
           y = prob, 
           col = variable, 
           shape = variable, 
           linetype = variable)) + 
  geom_line(size = 1.1) + 
  xlab(" ") + 
  ylab("%") + 
  geom_point(size = 2.5) +
  theme(panel.grid.minor = element_blank( ), 
        panel.grid.major = element_blank( ), 
        panel.background = element_rect(fill="white", 
                                        colour="grey50"), 
        legend.title = element_blank( ),
        legend.position = "bottom") +
  scale_x_continuous(breaks = seq(2009, 2018 , 1), 
                     labels = c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018")) + 
  scale_y_continuous(breaks = seq(0.00, 0.30, 0.05)) +
  expand_limits(y=c(0.00, 0.30)) #+
#labs(title = "Prevalences of C, E, and D in our simulated dataset at each T = t")

setwd("/Users/luissegura/Dropbox/Causal Trends/Figures")

ggsave("scenario2_fig1.jpg", width=9 ,height=6, units= "in", dpi=1200)


#--------------------------------------------------------#
#                                                        #
#   MODELING THE EFFECT OF C on E and D, and E on D      #
#                                                        #
#--------------------------------------------------------#

### FUNCTION TO SELECT A DF FROM THE LIST OF 100 DFs
list_select<-function(some_list, x){
  return(some_list[[x]])
}

### FUNCTION to extract coefficients
list_select_coeff<-function(some_list, x){
  return(some_list[[x]][["coefficients"]][[2]])
}

### Generating an object with all 100 datasets appended
df.all<-data.frame(matrix(NA, 1, 1)) ### generate an empty dataframe to store all the datasets appended

df.all<-list_select(listofdf, 1) ### select the first dataframe from the list and store it in the object

for(i in 2:100){
  df.all<-rbind(df.all, list_select(listofdf, i)) ### loop the selection of elements of a list and append them into the dataframe
}


df.all[, 7]<-rep(1:100, each = 10000)
colnames(df.all)[7]<-"simulation"

lista<-list()
listb<-list()
listc<-list()
listd<-list()

### Estimating the OR E | C, OR D |C and OR D | E for each dataset 1:100
for(i in 1:100){
  a<-glm(as.factor(exposure)~ as.factor(confounder), data = subset(df.all, simulation==i), family = quasibinomial(link = "logit"))
  lista[[i]]<-a
  b<-glm(as.factor(disease)~ as.factor(confounder), data = subset(df.all, simulation==i), family = quasibinomial(link = "logit"))
  listb[[i]]<-b
  c<-glm(as.factor(disease)~ as.factor(exposure), data = subset(df.all, simulation==i), family = quasibinomial(link = "logit"))
  listc[[i]]<-c
  d<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, simulation==i), 
         family = quasibinomial(link = "logit"))
  listd[[i]]<-d
}

### Extracting the coefficient OR E | C
for(i in 1:100){ 
  e_c_or<-exp(list_select_coeff(lista, i))
}

### Extracting the 95%CI OR E | C
e_c_ci<-data.frame(matrix(NA, nrow=1, ncol=2))
e_c_ci<-exp(confint(list_select(lista, 1))[2, ])
for(i in 2:100){
  e_c_ci<-data.frame(rbind(e_c_ci, exp(confint(list_select(lista, i))[2, ])))
}
rownames(e_c_ci)<-1:100
colnames(e_c_ci)<-c("lbci", "ubci")

### Extracting the coefficient OR D | C
for(i in 1:100){ 
  d_c_or<-exp(list_select_coeff(listb, i))
}

### Extracting the 95%CI OR D | C
d_c_ci<-data.frame(matrix(NA, nrow=1, ncol=2))
d_c_ci<-exp(confint(list_select(listb, 1))[2, ])
for(i in 2:100){
  d_c_ci<-data.frame(rbind(d_c_ci, exp(confint(list_select(listb, i))[2, ])))
}
rownames(d_c_ci)<-1:100
colnames(d_c_ci)<-c("lbci", "ubci")

### Extracting the coefficient OR D | E
for(i in 1:100){ 
  d_e_or<-exp(list_select_coeff(listc, i))
}

### Extracting the 95%CI OR D | E
d_e_ci<-data.frame(matrix(NA, nrow=1, ncol=2))
d_e_ci<-exp(confint(list_select(listc, 1))[2, ])
for(i in 2:100){
  d_e_ci<-data.frame(rbind(d_e_ci, exp(confint(list_select(listc, i))[2, ])))
}
rownames(d_e_ci)<-1:100
colnames(d_e_ci)<-c("lbci", "ubci")

### Extracting the coefficient OR D | E
for(i in 1:100){ 
  d_e_c_or<-exp(list_select_coeff(listd, i))
}

### Extracting the 95%CI OR D | E
d_e_c_ci<-data.frame(matrix(NA, nrow=1, ncol=2))
d_e_c_ci<-exp(confint(list_select(listd, 1))[2, ])
for(i in 2:100){
  d_e_c_ci<-data.frame(rbind(d_e_c_ci, exp(confint(list_select(listd, i))[2, ])))
}
rownames(d_e_c_ci)<-1:100
colnames(d_e_c_ci)<-c("lbci", "ubci")


df.overall<-data.frame(cbind(mean(e_c_or), mean(e_c_ci$lbci), mean(e_c_ci$ubci), 
                             mean(d_c_or), mean(d_c_ci$lbci), mean(d_c_ci$ubci), 
                             mean(d_e_or), mean(d_e_ci$lbci), mean(d_e_ci$ubci), 
                             mean(d_e_c_or), mean(d_e_c_ci$lbci), mean(d_e_c_ci$ubci)))

colnames(df.overall)<-c("OR_E_C", "LBCI", "UBCI", "OR_D_C", "LBCI", "UBCI", "OR_D_E", "LBCI", "UBCI", "OR_D_E_C", "LBCI", "UBCI")
pander(df.overall, caption = "Overall effect of C on E and D, and E on D")


#-------------------------------------------------------------------------------#
#                                                                               #
#        Estimate the OR E | C for each year from 100 simulated datasets        #
#                                                                               #
#-------------------------------------------------------------------------------#


### List of models for 2009
listofm09<-list()

for(i in 1:100){
  x<-glm(exposure~ confounder, data = subset(df.all, year==2009 & simulation==i), family = quasibinomial(link = "logit"))
  listofm09[[i]]<-x
}

### List of models for 2010
listofm10<-list()

for(i in 1:100){
  x<-glm(exposure~ confounder, data = subset(df.all, year==2010 & simulation==i), family = quasibinomial(link = "logit"))
  listofm10[[i]]<-x
}

### List of models for 2011
listofm11<-list()
for(i in 1:100){
  x<-glm(exposure~ confounder, data = subset(df.all, year==2011 & simulation==i), family = quasibinomial(link = "logit"))
  listofm11[[i]]<-x
}

### List of models for 2012
listofm12<-list()
for(i in 1:100){
  x<-glm(exposure~ confounder, data = subset(df.all, year==2012 & simulation==i), family = quasibinomial(link = "logit"))
  listofm12[[i]]<-x
}

### List of models for 2013
listofm13<-list()
for(i in 1:100){
  x<-glm(exposure~ confounder, data = subset(df.all, year==2013 & simulation==i), family = quasibinomial(link = "logit"))
  listofm13[[i]]<-x
}

### List of models for 2014
listofm14<-list()
for(i in 1:100){
  x<-glm(exposure~ confounder, data = subset(df.all, year==2014 & simulation==i), family = quasibinomial(link = "logit"))
  listofm14[[i]]<-x
}

### List of models for 2015
listofm15<-list()
for(i in 1:100){
  x<-glm(exposure~ confounder, data = subset(df.all, year==2015 & simulation==i), family = quasibinomial(link = "logit"))
  listofm15[[i]]<-x
}

### List of models for 2016
listofm16<-list()
for(i in 1:100){
  x<-glm(exposure~ confounder, data = subset(df.all, year==2016 & simulation==i), family = quasibinomial(link = "logit"))
  listofm16[[i]]<-x
}

### List of models for 2017
listofm17<-list()
for(i in 1:100){
  x<-glm(exposure~ confounder, data = subset(df.all, year==2017 & simulation==i), family = quasibinomial(link = "logit"))
  listofm17[[i]]<-x
}

### List of models for 2018
listofm18<-list()
for(i in 1:100){
  x<-glm(exposure~ confounder, data = subset(df.all, year==2018 & simulation==i), family = quasibinomial(link = "logit"))
  listofm18[[i]]<-x
}



### Extract coefficients for m09
or_09<-data.frame(matrix(NA, 1, 1))
or_09<-exp(listofm09[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2009 for the 100 simulations
for(i in 2:100){
  or_09[i]<-exp(list_select_coeff(listofm09, i))
}

### Extract coefficients for m10
or_10<-data.frame(matrix(NA, 1, 1))
or_10<-exp(listofm10[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2010 for the 100 simulations
for(i in 2:100){
  or_10[i]<-exp(list_select_coeff(listofm10, i))
}

### Extract coefficients for m11
or_11<-data.frame(matrix(NA, 1, 1))
or_11<-exp(listofm11[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2011 for the 100 simulations
for(i in 2:100){
  or_11[i]<-exp(list_select_coeff(listofm11, i))
}

### Extract coefficients for m12
or_12<-data.frame(matrix(NA, 1, 1))
or_12<-exp(listofm12[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2012 for the 100 simulations
for(i in 2:100){
  or_12[i]<-exp(list_select_coeff(listofm12, i))
}

### Extract coefficients for m13
or_13<-data.frame(matrix(NA, 1, 1))
or_13<-exp(listofm13[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2013 for the 100 simulations
for(i in 2:100){
  or_13[i]<-exp(list_select_coeff(listofm13, i))
}

### Extract coefficients for m14
or_14<-data.frame(matrix(NA, 1, 1))
or_14<-exp(listofm14[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2014 for the 100 simulations
for(i in 2:100){
  or_14[i]<-exp(list_select_coeff(listofm14, i))
}

### Extract coefficients for m15
or_15<-data.frame(matrix(NA, 1, 1))
or_15<-exp(listofm15[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2015 for the 100 simulations
for(i in 2:100){
  or_15[i]<-exp(list_select_coeff(listofm15, i))
}

### Extract coefficients for m16
or_16<-data.frame(matrix(NA, 1, 1))
or_16<-exp(listofm16[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2016 for the 100 simulations
for(i in 2:100){
  or_16[i]<-exp(list_select_coeff(listofm16, i))
}

### Extract coefficients for m17
or_17<-data.frame(matrix(NA, 1, 1))
or_17<-exp(listofm17[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2017 for the 100 simulations
for(i in 2:100){
  or_17[i]<-exp(list_select_coeff(listofm17, i))
}

### Extract coefficients for m18
or_18<-data.frame(matrix(NA, 1, 1))
or_18<-exp(listofm18[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2018 for the 100 simulations
for(i in 2:100){
  or_18[i]<-exp(list_select_coeff(listofm18, i))
}

### Extracting confint of OR for 2009
df.ci.09<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.09<-exp(confint(list_select(listofm09, 1))[2, ])

for(i in 2:100){
  df.ci.09<-data.frame(rbind(df.ci.09, exp(confint(list_select(listofm09, i))[2, ])))
}
rownames(df.ci.09)<-1:100
colnames(df.ci.09)<-c("lbci", "ubci")


### Extracting confint of OR for 2010
df.ci.10<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.10<-exp(confint(list_select(listofm10, 1))[2, ])

for(i in 2:100){
  df.ci.10<-data.frame(rbind(df.ci.10, exp(confint(list_select(listofm10, i))[2, ])))
}
rownames(df.ci.10)<-1:100
colnames(df.ci.10)<-c("lbci", "ubci")

### Extracting confint of OR for 2011
df.ci.11<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.11<-exp(confint(list_select(listofm11, 1))[2, ])

for(i in 2:100){
  df.ci.11<-data.frame(rbind(df.ci.11, exp(confint(list_select(listofm11, i))[2, ])))
}
rownames(df.ci.11)<-1:100
colnames(df.ci.11)<-c("lbci", "ubci")

### Extracting confint of OR for 2012
df.ci.12<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.12<-exp(confint(list_select(listofm12, 1))[2, ])

for(i in 2:100){
  df.ci.12<-data.frame(rbind(df.ci.12, exp(confint(list_select(listofm12, i))[2, ])))
}
rownames(df.ci.12)<-1:100
colnames(df.ci.12)<-c("lbci", "ubci")


### Extracting confint of OR for 2013
df.ci.13<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.13<-exp(confint(list_select(listofm13, 1))[2, ])

for(i in 2:100){
  df.ci.13<-data.frame(rbind(df.ci.13, exp(confint(list_select(listofm13, i))[2, ])))
}
rownames(df.ci.13)<-1:100
colnames(df.ci.13)<-c("lbci", "ubci")


### Extracting confint of OR for 2014
df.ci.14<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.14<-exp(confint(list_select(listofm14, 1))[2, ])

for(i in 2:100){
  df.ci.14<-data.frame(rbind(df.ci.14, exp(confint(list_select(listofm14, i))[2, ])))
}
rownames(df.ci.14)<-1:100
colnames(df.ci.14)<-c("lbci", "ubci")

### Extracting confint of OR for 2015
df.ci.15<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.15<-exp(confint(list_select(listofm15, 1))[2, ])

for(i in 2:100){
  df.ci.15<-data.frame(rbind(df.ci.15, exp(confint(list_select(listofm15, i))[2, ])))
}
rownames(df.ci.15)<-1:100
colnames(df.ci.15)<-c("lbci", "ubci")

### Extracting confint of OR for 2016
df.ci.16<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.16<-exp(confint(list_select(listofm16, 1))[2, ])

for(i in 2:100){
  df.ci.16<-data.frame(rbind(df.ci.16, exp(confint(list_select(listofm16, i))[2, ])))
}
rownames(df.ci.16)<-1:100
colnames(df.ci.16)<-c("lbci", "ubci")

### Extracting confint of OR for 2017
df.ci.17<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.17<-exp(confint(list_select(listofm17, 1))[2, ])

for(i in 2:100){
  df.ci.17<-data.frame(rbind(df.ci.17, exp(confint(list_select(listofm17, i))[2, ])))
}
rownames(df.ci.17)<-1:100
colnames(df.ci.17)<-c("lbci", "ubci")

### Extracting confint of OR for 2018
df.ci.18<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.18<-exp(confint(list_select(listofm18, 1))[2, ])

for(i in 2:100){
  df.ci.18<-data.frame(rbind(df.ci.18, exp(confint(list_select(listofm18, i))[2, ])))
}
rownames(df.ci.18)<-1:100
colnames(df.ci.18)<-c("lbci", "ubci")


df.all.coef1<-data.frame(year=2009:2018, 
                         or_e_c= rbind(mean(or_09), mean(or_10), mean(or_11), mean(or_12), mean(or_13), 
                                       mean(or_14), mean(or_15), mean(or_16), mean(or_17), mean(or_18)), 
                         lbci_e_c= rbind(mean(df.ci.09$lbci), mean(df.ci.10$lbci), mean(df.ci.11$lbci), mean(df.ci.12$lbci), 
                                         mean(df.ci.13$lbci), mean(df.ci.14$lbci), mean(df.ci.15$lbci), mean(df.ci.16$lbci), 
                                         mean(df.ci.17$lbci), mean(df.ci.18$lbci)), 
                         ubci_e_c= rbind(mean(df.ci.09$ubci), mean(df.ci.10$ubci), mean(df.ci.11$ubci), mean(df.ci.12$ubci), 
                                         mean(df.ci.13$ubci), mean(df.ci.14$ubci), mean(df.ci.15$ubci), mean(df.ci.16$ubci), 
                                         mean(df.ci.17$ubci), mean(df.ci.18$ubci)))

pander(df.all.coef1, caption = "Change over time in the effect (OR) of E | C using 100 simulated datasets")


#-------------------------------------------------------------------------------#
#                                                                               #
#        Estimate the OR D | C for each year from 100 simulated datasets        #
#                                                                               #
#-------------------------------------------------------------------------------#

### List of models for 2009
listofm209<-list()

for(i in 1:100){
  x<-glm(disease~ confounder, data = subset(df.all, year==2009 & simulation==i), family = quasibinomial(link = "logit"))
  listofm209[[i]]<-x
}

### List of models for 2010
listofm210<-list()

for(i in 1:100){
  x<-glm(disease~ confounder, data = subset(df.all, year==2010 & simulation==i), family = quasibinomial(link = "logit"))
  listofm210[[i]]<-x
}

### List of models for 2011
listofm211<-list()
for(i in 1:100){
  x<-glm(disease~ confounder, data = subset(df.all, year==2011 & simulation==i), family = quasibinomial(link = "logit"))
  listofm211[[i]]<-x
}

### List of models for 2012
listofm212<-list()
for(i in 1:100){
  x<-glm(disease~ confounder, data = subset(df.all, year==2012 & simulation==i), family = quasibinomial(link = "logit"))
  listofm212[[i]]<-x
}

### List of models for 2013
listofm213<-list()
for(i in 1:100){
  x<-glm(disease~ confounder, data = subset(df.all, year==2013 & simulation==i), family = quasibinomial(link = "logit"))
  listofm213[[i]]<-x
}

### List of models for 2014
listofm214<-list()
for(i in 1:100){
  x<-glm(disease~ confounder, data = subset(df.all, year==2014 & simulation==i), family = quasibinomial(link = "logit"))
  listofm214[[i]]<-x
}

### List of models for 2015
listofm215<-list()
for(i in 1:100){
  x<-glm(disease~ confounder, data = subset(df.all, year==2015 & simulation==i), family = quasibinomial(link = "logit"))
  listofm215[[i]]<-x
}

### List of models for 2016
listofm216<-list()
for(i in 1:100){
  x<-glm(disease~ confounder, data = subset(df.all, year==2016 & simulation==i), family = quasibinomial(link = "logit"))
  listofm216[[i]]<-x
}

### List of models for 2017
listofm217<-list()
for(i in 1:100){
  x<-glm(disease~ confounder, data = subset(df.all, year==2017 & simulation==i), family = quasibinomial(link = "logit"))
  listofm217[[i]]<-x
}

### List of models for 2018
listofm218<-list()
for(i in 1:100){
  x<-glm(disease~ confounder, data = subset(df.all, year==2018 & simulation==i), family = quasibinomial(link = "logit"))
  listofm218[[i]]<-x
}

### Extract coefficients for m09
or_209<-data.frame(matrix(NA, 1, 1))
or_209<-exp(listofm209[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2009 for the 100 simulations
for(i in 2:100){
  or_209[i]<-exp(list_select_coeff(listofm209, i))
}

### Extract coefficients for m10
or_210<-data.frame(matrix(NA, 1, 1))
or_210<-exp(listofm210[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2010 for the 100 simulations
for(i in 2:100){
  or_210[i]<-exp(list_select_coeff(listofm210, i))
}

### Extract coefficients for m11
or_211<-data.frame(matrix(NA, 1, 1))
or_211<-exp(listofm211[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2011 for the 100 simulations
for(i in 2:100){
  or_211[i]<-exp(list_select_coeff(listofm211, i))
}

### Extract coefficients for m12
or_212<-data.frame(matrix(NA, 1, 1))
or_212<-exp(listofm212[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2012 for the 100 simulations
for(i in 2:100){
  or_212[i]<-exp(list_select_coeff(listofm212, i))
}

### Extract coefficients for m13
or_213<-data.frame(matrix(NA, 1, 1))
or_213<-exp(listofm213[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2013 for the 100 simulations
for(i in 2:100){
  or_213[i]<-exp(list_select_coeff(listofm213, i))
}

### Extract coefficients for m14
or_214<-data.frame(matrix(NA, 1, 1))
or_214<-exp(listofm214[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2014 for the 100 simulations
for(i in 2:100){
  or_214[i]<-exp(list_select_coeff(listofm214, i))
}

### Extract coefficients for m15
or_215<-data.frame(matrix(NA, 1, 1))
or_215<-exp(listofm215[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2015 for the 100 simulations
for(i in 2:100){
  or_215[i]<-exp(list_select_coeff(listofm215, i))
}

### Extract coefficients for m16
or_216<-data.frame(matrix(NA, 1, 1))
or_216<-exp(listofm216[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2016 for the 100 simulations
for(i in 2:100){
  or_216[i]<-exp(list_select_coeff(listofm216, i))
}

### Extract coefficients for m17
or_217<-data.frame(matrix(NA, 1, 1))
or_217<-exp(listofm217[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2017 for the 100 simulations
for(i in 2:100){
  or_217[i]<-exp(list_select_coeff(listofm217, i))
}

### Extract coefficients for m18
or_218<-data.frame(matrix(NA, 1, 1))
or_218<-exp(listofm218[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2018 for the 100 simulations
for(i in 2:100){
  or_218[i]<-exp(list_select_coeff(listofm218, i))
}

### Extracting confint of OR for 2009
df.ci.209<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.209<-exp(confint(list_select(listofm209, 1))[2, ])

for(i in 2:100){
  df.ci.209<-data.frame(rbind(df.ci.209, exp(confint(list_select(listofm209, i))[2, ])))
}
rownames(df.ci.209)<-1:100
colnames(df.ci.209)<-c("lbci", "ubci")


### Extracting confint of OR for 2010
df.ci.210<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.210<-exp(confint(list_select(listofm210, 1))[2, ])

for(i in 2:100){
  df.ci.210<-data.frame(rbind(df.ci.210, exp(confint(list_select(listofm210, i))[2, ])))
}
rownames(df.ci.210)<-1:100
colnames(df.ci.210)<-c("lbci", "ubci")

### Extracting confint of OR for 2011
df.ci.211<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.211<-exp(confint(list_select(listofm211, 1))[2, ])

for(i in 2:100){
  df.ci.211<-data.frame(rbind(df.ci.211, exp(confint(list_select(listofm211, i))[2, ])))
}
rownames(df.ci.211)<-1:100
colnames(df.ci.211)<-c("lbci", "ubci")

### Extracting confint of OR for 2012
df.ci.212<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.212<-exp(confint(list_select(listofm212, 1))[2, ])

for(i in 2:100){
  df.ci.212<-data.frame(rbind(df.ci.212, exp(confint(list_select(listofm212, i))[2, ])))
}
rownames(df.ci.212)<-1:100
colnames(df.ci.212)<-c("lbci", "ubci")


### Extracting confint of OR for 2013
df.ci.213<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.213<-exp(confint(list_select(listofm213, 1))[2, ])

for(i in 2:100){
  df.ci.213<-data.frame(rbind(df.ci.213, exp(confint(list_select(listofm213, i))[2, ])))
}
rownames(df.ci.213)<-1:100
colnames(df.ci.213)<-c("lbci", "ubci")


### Extracting confint of OR for 2014
df.ci.214<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.214<-exp(confint(list_select(listofm214, 1))[2, ])

for(i in 2:100){
  df.ci.214<-data.frame(rbind(df.ci.214, exp(confint(list_select(listofm214, i))[2, ])))
}
rownames(df.ci.214)<-1:100
colnames(df.ci.214)<-c("lbci", "ubci")

### Extracting confint of OR for 2015
df.ci.215<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.215<-exp(confint(list_select(listofm215, 1))[2, ])

for(i in 2:100){
  df.ci.215<-data.frame(rbind(df.ci.215, exp(confint(list_select(listofm215, i))[2, ])))
}
rownames(df.ci.215)<-1:100
colnames(df.ci.215)<-c("lbci", "ubci")

### Extracting confint of OR for 2016
df.ci.216<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.216<-exp(confint(list_select(listofm216, 1))[2, ])

for(i in 2:100){
  df.ci.216<-data.frame(rbind(df.ci.216, exp(confint(list_select(listofm216, i))[2, ])))
}
rownames(df.ci.216)<-1:100
colnames(df.ci.216)<-c("lbci", "ubci")

### Extracting confint of OR for 2017
df.ci.217<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.217<-exp(confint(list_select(listofm217, 1))[2, ])

for(i in 2:100){
  df.ci.217<-data.frame(rbind(df.ci.217, exp(confint(list_select(listofm217, i))[2, ])))
}
rownames(df.ci.217)<-1:100
colnames(df.ci.217)<-c("lbci", "ubci")

### Extracting confint of OR for 2018
df.ci.218<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.218<-exp(confint(list_select(listofm218, 1))[2, ])

for(i in 2:100){
  df.ci.218<-data.frame(rbind(df.ci.218, exp(confint(list_select(listofm218, i))[2, ])))
}
rownames(df.ci.218)<-1:100
colnames(df.ci.218)<-c("lbci", "ubci")


df.all.coef2<-data.frame(year=2009:2018, 
                         or_d_c= rbind(mean(or_209), mean(or_210), mean(or_211), mean(or_212), mean(or_213), 
                                       mean(or_214), mean(or_215), mean(or_216), mean(or_217), mean(or_218)), 
                         lbci_d_c= rbind(mean(df.ci.209$lbci), mean(df.ci.210$lbci), mean(df.ci.211$lbci), mean(df.ci.212$lbci), 
                                         mean(df.ci.213$lbci), mean(df.ci.214$lbci), mean(df.ci.215$lbci), mean(df.ci.216$lbci), 
                                         mean(df.ci.217$lbci), mean(df.ci.218$lbci)), 
                         ubci_d_c= rbind(mean(df.ci.209$ubci), mean(df.ci.210$ubci), mean(df.ci.211$ubci), mean(df.ci.212$ubci), 
                                         mean(df.ci.213$ubci), mean(df.ci.214$ubci), mean(df.ci.215$ubci), mean(df.ci.216$ubci), 
                                         mean(df.ci.217$ubci), mean(df.ci.218$ubci)))

pander(df.all.coef2, caption = "Change over time in the effect (OR) of D | C using 100 simulated datasets")


#-------------------------------------------------------------------------------#
#                                                                               #
#        Estimate the OR D | E for each year from 100 simulated datasets        #
#                                                                               #
#-------------------------------------------------------------------------------#

### List of models for 2009
listofm309<-list()

for(i in 1:100){
  x<-glm(disease~ exposure, data = subset(df.all, year==2009 & simulation==i), family = quasibinomial(link = "logit"))
  listofm309[[i]]<-x
}

### List of models for 2010
listofm310<-list()

for(i in 1:100){
  x<-glm(disease~ exposure, data = subset(df.all, year==2010 & simulation==i), family = quasibinomial(link = "logit"))
  listofm310[[i]]<-x
}

### List of models for 2011
listofm311<-list()

for(i in 1:100){
  x<-glm(disease~ exposure, data = subset(df.all, year==2011 & simulation==i), family = quasibinomial(link = "logit"))
  listofm311[[i]]<-x
}

### List of models for 2012
listofm312<-list()

for(i in 1:100){
  x<-glm(disease~ exposure, data = subset(df.all, year==2012 & simulation==i), family = quasibinomial(link = "logit"))
  listofm312[[i]]<-x
}

### List of models for 2013
listofm313<-list()

for(i in 1:100){
  x<-glm(disease~ exposure, data = subset(df.all, year==2013 & simulation==i), family = quasibinomial(link = "logit"))
  listofm313[[i]]<-x
}

### List of models for 2014
listofm314<-list()

for(i in 1:100){
  x<-glm(disease~ exposure, data = subset(df.all, year==2014 & simulation==i), family = quasibinomial(link = "logit"))
  listofm314[[i]]<-x
}

### List of models for 2015
listofm315<-list()
for(i in 1:100){
  x<-glm(disease~ exposure, data = subset(df.all, year==2015 & simulation==i), family = quasibinomial(link = "logit"))
  listofm315[[i]]<-x
}

### List of models for 2016
listofm316<-list()
for(i in 1:100){
  x<-glm(disease~ exposure, data = subset(df.all, year==2016 & simulation==i), family = quasibinomial(link = "logit"))
  listofm316[[i]]<-x
}

### List of models for 2017
listofm317<-list()
for(i in 1:100){
  x<-glm(disease~ exposure, data = subset(df.all, year==2017 & simulation==i), family = quasibinomial(link = "logit"))
  listofm317[[i]]<-x
}

### List of models for 2018
listofm318<-list()
for(i in 1:100){
  x<-glm(disease~ exposure, data = subset(df.all, year==2018 & simulation==i), family = quasibinomial(link = "logit"))
  listofm318[[i]]<-x
}

### Extract coefficients for m09
or_309<-data.frame(matrix(NA, 1, 1))
or_309<-exp(listofm309[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2009 for the 100 simulations
for(i in 2:100){
  or_309[i]<-exp(list_select_coeff(listofm309, i))
}

### Extract coefficients for m10
or_310<-data.frame(matrix(NA, 1, 1))
or_310<-exp(listofm310[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2010 for the 100 simulations
for(i in 2:100){
  or_310[i]<-exp(list_select_coeff(listofm310, i))
}

### Extract coefficients for m11
or_311<-data.frame(matrix(NA, 1, 1))
or_311<-exp(listofm311[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2011 for the 100 simulations
for(i in 2:100){
  or_311[i]<-exp(list_select_coeff(listofm311, i))
}

### Extract coefficients for m12
or_312<-data.frame(matrix(NA, 1, 1))
or_312<-exp(listofm312[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2012 for the 100 simulations
for(i in 2:100){
  or_312[i]<-exp(list_select_coeff(listofm312, i))
}

### Extract coefficients for m13
or_313<-data.frame(matrix(NA, 1, 1))
or_313<-exp(listofm313[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2013 for the 100 simulations
for(i in 2:100){
  or_313[i]<-exp(list_select_coeff(listofm313, i))
}

### Extract coefficients for m14
or_314<-data.frame(matrix(NA, 1, 1))
or_314<-exp(listofm314[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2014 for the 100 simulations
for(i in 2:100){
  or_314[i]<-exp(list_select_coeff(listofm314, i))
}

### Extract coefficients for m15
or_315<-data.frame(matrix(NA, 1, 1))
or_315<-exp(listofm315[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2015 for the 100 simulations
for(i in 2:100){
  or_315[i]<-exp(list_select_coeff(listofm315, i))
}

### Extract coefficients for m16
or_316<-data.frame(matrix(NA, 1, 1))
or_316<-exp(listofm316[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2016 for the 100 simulations
for(i in 2:100){
  or_316[i]<-exp(list_select_coeff(listofm316, i))
}

### Extract coefficients for m17
or_317<-data.frame(matrix(NA, 1, 1))
or_317<-exp(listofm317[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2017 for the 100 simulations
for(i in 2:100){
  or_317[i]<-exp(list_select_coeff(listofm317, i))
}

### Extract coefficients for m18
or_318<-data.frame(matrix(NA, 1, 1))
or_318<-exp(listofm318[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2018 for the 100 simulations
for(i in 2:100){
  or_318[i]<-exp(list_select_coeff(listofm318, i))
}


### Extracting confint of OR for 2009
df.ci.309<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.309<-exp(confint(list_select(listofm309, 1))[2, ])

for(i in 2:100){
  df.ci.309<-data.frame(rbind(df.ci.309, exp(confint(list_select(listofm309, i))[2, ])))
}
rownames(df.ci.309)<-1:100
colnames(df.ci.309)<-c("lbci", "ubci")


### Extracting confint of OR for 2010
df.ci.310<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.310<-exp(confint(list_select(listofm310, 1))[2, ])

for(i in 2:100){
  df.ci.310<-data.frame(rbind(df.ci.310, exp(confint(list_select(listofm310, i))[2, ])))
}
rownames(df.ci.310)<-1:100
colnames(df.ci.310)<-c("lbci", "ubci")

### Extracting confint of OR for 2011
df.ci.311<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.311<-exp(confint(list_select(listofm311, 1))[2, ])

for(i in 2:100){
  df.ci.311<-data.frame(rbind(df.ci.311, exp(confint(list_select(listofm311, i))[2, ])))
}
rownames(df.ci.311)<-1:100
colnames(df.ci.311)<-c("lbci", "ubci")

### Extracting confint of OR for 2012
df.ci.312<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.312<-exp(confint(list_select(listofm312, 1))[2, ])

for(i in 2:100){
  df.ci.312<-data.frame(rbind(df.ci.312, exp(confint(list_select(listofm312, i))[2, ])))
}
rownames(df.ci.312)<-1:100
colnames(df.ci.312)<-c("lbci", "ubci")


### Extracting confint of OR for 2013
df.ci.313<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.313<-exp(confint(list_select(listofm313, 1))[2, ])

for(i in 2:100){
  df.ci.313<-data.frame(rbind(df.ci.313, exp(confint(list_select(listofm313, i))[2, ])))
}
rownames(df.ci.313)<-1:100
colnames(df.ci.313)<-c("lbci", "ubci")


### Extracting confint of OR for 2014
df.ci.314<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.314<-exp(confint(list_select(listofm314, 1))[2, ])

for(i in 2:100){
  df.ci.314<-data.frame(rbind(df.ci.314, exp(confint(list_select(listofm314, i))[2, ])))
}
rownames(df.ci.314)<-1:100
colnames(df.ci.314)<-c("lbci", "ubci")

### Extracting confint of OR for 2015
df.ci.315<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.315<-exp(confint(list_select(listofm315, 1))[2, ])

for(i in 2:100){
  df.ci.315<-data.frame(rbind(df.ci.315, exp(confint(list_select(listofm315, i))[2, ])))
}
rownames(df.ci.315)<-1:100
colnames(df.ci.315)<-c("lbci", "ubci")

### Extracting confint of OR for 2016
df.ci.316<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.316<-exp(confint(list_select(listofm316, 1))[2, ])

for(i in 2:100){
  df.ci.316<-data.frame(rbind(df.ci.316, exp(confint(list_select(listofm316, i))[2, ])))
}
rownames(df.ci.316)<-1:100
colnames(df.ci.316)<-c("lbci", "ubci")

### Extracting confint of OR for 2017
df.ci.317<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.317<-exp(confint(list_select(listofm317, 1))[2, ])

for(i in 2:100){
  df.ci.317<-data.frame(rbind(df.ci.317, exp(confint(list_select(listofm317, i))[2, ])))
}
rownames(df.ci.317)<-1:100
colnames(df.ci.317)<-c("lbci", "ubci")

### Extracting confint of OR for 2018
df.ci.318<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.318<-exp(confint(list_select(listofm318, 1))[2, ])

for(i in 2:100){
  df.ci.318<-data.frame(rbind(df.ci.318, exp(confint(list_select(listofm318, i))[2, ])))
}
rownames(df.ci.318)<-1:100
colnames(df.ci.318)<-c("lbci", "ubci")


df.all.coef3<-data.frame(year=2009:2018, 
                         or_d_e= rbind(mean(or_309), mean(or_310), mean(or_311), mean(or_312), mean(or_313), 
                                       mean(or_314), mean(or_315), mean(or_316), mean(or_317), mean(or_318)), 
                         lbci_d_e= rbind(mean(df.ci.309$lbci), mean(df.ci.310$lbci), mean(df.ci.311$lbci), mean(df.ci.312$lbci), 
                                         mean(df.ci.313$lbci), mean(df.ci.314$lbci), mean(df.ci.315$lbci), mean(df.ci.316$lbci), 
                                         mean(df.ci.317$lbci), mean(df.ci.318$lbci)), 
                         ubci_d_e= rbind(mean(df.ci.309$ubci), mean(df.ci.310$ubci), mean(df.ci.311$ubci), mean(df.ci.312$ubci), 
                                         mean(df.ci.313$ubci), mean(df.ci.314$ubci), mean(df.ci.315$ubci), mean(df.ci.316$ubci), 
                                         mean(df.ci.317$ubci), mean(df.ci.318$ubci)))

pander(df.all.coef3, caption = "Change over time in the effect (OR) of E | C using 100 simulated datasets")


#-------------------------------------------------------------------------------#
#                                                                               #
#     Estimate the OR D | E, C for each year from 100 simulated datasets        #
#                                                                               #
#-------------------------------------------------------------------------------#

### List of models for 2009
listofm409<-list()

for(i in 1:100){
  x<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, year==2009 & simulation==i), 
         family = quasibinomial(link = "logit"))
  listofm409[[i]]<-x
}

### List of models for 2010
listofm410<-list()

for(i in 1:100){
  x<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, year==2010 & simulation==i), 
         family = quasibinomial(link = "logit"))
  listofm410[[i]]<-x
}

### List of models for 2011
listofm411<-list()

for(i in 1:100){
  x<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, year==2011 & simulation==i), 
         family = quasibinomial(link = "logit"))
  listofm411[[i]]<-x
}

### List of models for 2012
listofm412<-list()

for(i in 1:100){
  x<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, year==2012 & simulation==i), 
         family = quasibinomial(link = "logit"))
  listofm412[[i]]<-x
}

### List of models for 2013
listofm413<-list()

for(i in 1:100){
  x<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, year==2013 & simulation==i), 
         family = quasibinomial(link = "logit"))
  listofm413[[i]]<-x
}

### List of models for 2014
listofm414<-list()

for(i in 1:100){
  x<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, year==2014 & simulation==i), 
         family = quasibinomial(link = "logit"))
  listofm414[[i]]<-x
}

### List of models for 2015
listofm415<-list()
for(i in 1:100){
  x<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, year==2015 & simulation==i), 
         family = quasibinomial(link = "logit"))
  listofm415[[i]]<-x
}

### List of models for 2016
listofm416<-list()
for(i in 1:100){
  x<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, year==2016 & simulation==i), 
         family = quasibinomial(link = "logit"))
  listofm416[[i]]<-x
}

### List of models for 2017
listofm417<-list()
for(i in 1:100){
  x<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, year==2017 & simulation==i), 
         family = quasibinomial(link = "logit"))
  listofm417[[i]]<-x
}

### List of models for 2018
listofm418<-list()
for(i in 1:100){
  x<-glm(as.factor(disease)~ as.factor(exposure) + as.factor(confounder), data = subset(df.all, year==2018 & simulation==i), 
         family = quasibinomial(link = "logit"))
  listofm418[[i]]<-x
}


### Extract coefficients for m09
or_409<-data.frame(matrix(NA, 1, 1))
or_409<-exp(listofm409[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2009 for the 100 simulations
for(i in 2:100){
  or_409[i]<-exp(list_select_coeff(listofm409, i))
}

### Extract coefficients for m10
or_410<-data.frame(matrix(NA, 1, 1))
or_410<-exp(listofm410[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2010 for the 100 simulations
for(i in 2:100){
  or_410[i]<-exp(list_select_coeff(listofm410, i))
}

### Extract coefficients for m11
or_411<-data.frame(matrix(NA, 1, 1))
or_411<-exp(listofm411[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2011 for the 100 simulations
for(i in 2:100){
  or_411[i]<-exp(list_select_coeff(listofm411, i))
}

### Extract coefficients for m12
or_412<-data.frame(matrix(NA, 1, 1))
or_412<-exp(listofm412[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2012 for the 100 simulations
for(i in 2:100){
  or_412[i]<-exp(list_select_coeff(listofm412, i))
}

### Extract coefficients for m13
or_413<-data.frame(matrix(NA, 1, 1))
or_413<-exp(listofm413[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2013 for the 100 simulations
for(i in 2:100){
  or_413[i]<-exp(list_select_coeff(listofm413, i))
}

### Extract coefficients for m14
or_414<-data.frame(matrix(NA, 1, 1))
or_414<-exp(listofm414[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2014 for the 100 simulations
for(i in 2:100){
  or_414[i]<-exp(list_select_coeff(listofm414, i))
}

### Extract coefficients for m15
or_415<-data.frame(matrix(NA, 1, 1))
or_415<-exp(listofm415[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2015 for the 100 simulations
for(i in 2:100){
  or_415[i]<-exp(list_select_coeff(listofm415, i))
}

### Extract coefficients for m16
or_416<-data.frame(matrix(NA, 1, 1))
or_416<-exp(listofm416[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2016 for the 100 simulations
for(i in 2:100){
  or_416[i]<-exp(list_select_coeff(listofm416, i))
}

### Extract coefficients for m17
or_417<-data.frame(matrix(NA, 1, 1))
or_417<-exp(listofm417[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2017 for the 100 simulations
for(i in 2:100){
  or_417[i]<-exp(list_select_coeff(listofm417, i))
}

### Extract coefficients for m18
or_418<-data.frame(matrix(NA, 1, 1))
or_418<-exp(listofm418[[1]][["coefficients"]][[2]])

### Loop it to create a dataframe with all ORs of 2018 for the 100 simulations
for(i in 2:100){
  or_418[i]<-exp(list_select_coeff(listofm418, i))
}

### Extracting confint of OR for 2009
df.ci.409<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.409<-exp(confint(list_select(listofm409, 1))[2, ])

for(i in 2:100){
  df.ci.409<-data.frame(rbind(df.ci.409, exp(confint(list_select(listofm409, i))[2, ])))
}
rownames(df.ci.409)<-1:100
colnames(df.ci.409)<-c("lbci", "ubci")


### Extracting confint of OR for 2010
df.ci.410<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.410<-exp(confint(list_select(listofm410, 1))[2, ])

for(i in 2:100){
  df.ci.410<-data.frame(rbind(df.ci.410, exp(confint(list_select(listofm410, i))[2, ])))
}
rownames(df.ci.410)<-1:100
colnames(df.ci.410)<-c("lbci", "ubci")

### Extracting confint of OR for 2011
df.ci.411<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.411<-exp(confint(list_select(listofm411, 1))[2, ])

for(i in 2:100){
  df.ci.411<-data.frame(rbind(df.ci.411, exp(confint(list_select(listofm411, i))[2, ])))
}
rownames(df.ci.411)<-1:100
colnames(df.ci.411)<-c("lbci", "ubci")

### Extracting confint of OR for 2012
df.ci.412<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.412<-exp(confint(list_select(listofm412, 1))[2, ])

for(i in 2:100){
  df.ci.412<-data.frame(rbind(df.ci.412, exp(confint(list_select(listofm412, i))[2, ])))
}
rownames(df.ci.412)<-1:100
colnames(df.ci.412)<-c("lbci", "ubci")


### Extracting confint of OR for 2013
df.ci.413<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.413<-exp(confint(list_select(listofm413, 1))[2, ])

for(i in 2:100){
  df.ci.413<-data.frame(rbind(df.ci.413, exp(confint(list_select(listofm413, i))[2, ])))
}
rownames(df.ci.413)<-1:100
colnames(df.ci.413)<-c("lbci", "ubci")


### Extracting confint of OR for 2014
df.ci.414<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.414<-exp(confint(list_select(listofm414, 1))[2, ])

for(i in 2:100){
  df.ci.414<-data.frame(rbind(df.ci.414, exp(confint(list_select(listofm414, i))[2, ])))
}
rownames(df.ci.414)<-1:100
colnames(df.ci.414)<-c("lbci", "ubci")

### Extracting confint of OR for 2015
df.ci.415<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.415<-exp(confint(list_select(listofm415, 1))[2, ])

for(i in 2:100){
  df.ci.415<-data.frame(rbind(df.ci.415, exp(confint(list_select(listofm415, i))[2, ])))
}
rownames(df.ci.415)<-1:100
colnames(df.ci.415)<-c("lbci", "ubci")

### Extracting confint of OR for 2016
df.ci.416<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.416<-exp(confint(list_select(listofm416, 1))[2, ])

for(i in 2:100){
  df.ci.416<-data.frame(rbind(df.ci.416, exp(confint(list_select(listofm416, i))[2, ])))
}
rownames(df.ci.416)<-1:100
colnames(df.ci.416)<-c("lbci", "ubci")

### Extracting confint of OR for 2017
df.ci.417<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.417<-exp(confint(list_select(listofm417, 1))[2, ])

for(i in 2:100){
  df.ci.417<-data.frame(rbind(df.ci.417, exp(confint(list_select(listofm417, i))[2, ])))
}
rownames(df.ci.417)<-1:100
colnames(df.ci.417)<-c("lbci", "ubci")

### Extracting confint of OR for 2018
df.ci.418<-data.frame(matrix(NA, nrow=1, ncol=2))
df.ci.418<-exp(confint(list_select(listofm418, 1))[2, ])

for(i in 2:100){
  df.ci.418<-data.frame(rbind(df.ci.418, exp(confint(list_select(listofm418, i))[2, ])))
}
rownames(df.ci.418)<-1:100
colnames(df.ci.418)<-c("lbci", "ubci")

df.all.coef4<-data.frame(year=2009:2018, 
                         or_d_e_c= rbind(mean(or_409), mean(or_410), mean(or_411), mean(or_412), mean(or_413), 
                                         mean(or_414), mean(or_415), mean(or_416), mean(or_417), mean(or_418)), 
                         lbci_d_e_c= rbind(mean(df.ci.409$lbci), mean(df.ci.410$lbci), mean(df.ci.411$lbci), mean(df.ci.412$lbci), 
                                           mean(df.ci.413$lbci), mean(df.ci.414$lbci), mean(df.ci.415$lbci), mean(df.ci.416$lbci), 
                                           mean(df.ci.417$lbci), mean(df.ci.418$lbci)), 
                         ubci_d_e_c= rbind(mean(df.ci.409$ubci), mean(df.ci.410$ubci), mean(df.ci.411$ubci), mean(df.ci.412$ubci), 
                                           mean(df.ci.413$ubci), mean(df.ci.414$ubci), mean(df.ci.415$ubci), mean(df.ci.416$ubci), 
                                           mean(df.ci.417$ubci), mean(df.ci.418$ubci)))

pander(df.all.coef4, caption = "Change over time in the effect (OR) of D | E, C using 100 simulated datasets")
df.all.coef4

df.all.coef<-cbind(df.all.coef1, df.all.coef2, df.all.coef3, df.all.coef4)
pander(df.all.coef, caption = "Changer over time in the effects of C on E and D, and E on D")

df.all.coef
df.plot3<-data.frame(matrix(NA, nrow=40, ncol=5))
df.plot3[, 1]<-df.all.coef$year
df.plot3[, 2]<-rep(c("OR E | C", "OR D | C", "OR D | E", "OR D | E, C"), each = 10)
df.plot3[ , 3]<-c(df.all.coef$or_e_c, df.all.coef$or_d_c, df.all.coef$or_d_e, df.all.coef$or_d_e_c)
df.plot3[ , 4]<-c(df.all.coef$lbci_e_c, df.all.coef$lbci_d_c, df.all.coef$lbci_d_e, df.all.coef$lbci_d_e_c)
df.plot3[ , 5]<-c(df.all.coef$ubci_e_c, df.all.coef$ubci_d_c, df.all.coef$ubci_d_e, df.all.coef$ubci_d_e_c)
colnames(df.plot3)<-c("year", "effect", "or","lbci", "ubci")


ggplot(df.plot3, aes(x=year, y=or, ymin=lbci, ymax=ubci, col=effect, shape=effect)) +
  geom_pointrange(size = .8) +
  geom_hline(yintercept = 1, lty = 2) +
  coord_flip() +
  xlab (" ") +
  ylab ("OR (95% CI)") + 
  #geom_point(size = 2.5) +
  theme(panel.grid.minor = element_blank( ), 
        panel.grid.major = element_blank( ), 
        panel.background = element_rect(fill="white", 
                                        colour="grey50"), 
        legend.title = element_blank( ),
        legend.position = "bottom") +
  scale_x_continuous(breaks = seq(2009, 2018 , 1), 
                     labels = c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018")) + 
  scale_y_continuous(breaks = seq(0.00, 28.00, 2)) +
  scale_shape_manual(values = c(15, 16, 17, 18))
expand_limits(y=c(0.00, 28.00)) #+
#labs(title = "Prevalences of C, E, and D in our simulated dataset at each T = t")

setwd("/Users/luissegura/Dropbox/Causal Trends/Figures")

ggsave("scenario2_fig2_vertical.jpg", width=9 ,height=6, units= "in", dpi=1200)

scemario2.or<-df.all.coef
