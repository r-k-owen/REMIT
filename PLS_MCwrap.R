rm(list=ls())
library(dplyr)
library(flexsurv)
library(dmutate)
#
#memory.limit(size=56000)

library(ggplot2)
library(tidyr)
library(survival)
library(scales)  # for percent_format

setwd("/Users/rkowen/Library/CloudStorage/OneDrive-SwanseaUniversity/Research_Projects/NIHR_REMIT/HTA_Sim/MonteCarlo/")

######################
##Simulate a cohort ##
######################
#age<-rnorm(650000,65,7)
age<-rnorm(1000000,60,7)
pop<-as.data.frame(age)
pop$agec<- pop$age-65
pop$id<- seq_len(nrow(pop))
pop <- pop[,c(3,1,2)]
pop$sex <- rbinom(nrow(pop), size=1, prob=0.5)
prop.table(pop$sex)

wimd_levels <- c(
  "LeastDep",
  "LowMidDep",
  "MiddleDep",
  "HighMidDep",
  "MostDep"
)

probs <- c(0.20, 0.20, 0.20, 0.20, 0.20)  # must sum to 1

pop$wimd <- factor(
  sample(wimd_levels, size = nrow(pop), replace = TRUE, prob = probs),
  levels = wimd_levels
)

table(pop$wimd)
prop.table(table(pop$wimd))

wimd_dummies <- sapply(
  wimd_levels,
  function(lvl) as.integer(pop$wimd == lvl)
)

colnames(wimd_dummies) <- wimd_levels

pop$wimd1<- wimd_dummies[,"MostDep"]
pop$wimd2<- wimd_dummies[,"HighMidDep"]
pop$wimd3<- wimd_dummies[,"MiddleDep"]
pop$wimd4<- wimd_dummies[,"LowMidDep"]
pop$wimd5<- wimd_dummies[,"LeastDep"]

set.seed(1234)

# Monte Carlo wrapper #
B <- 100   # number of simulations
n <- nrow(pop)

results <- vector("list", B)
mc_pop      <- vector("list", B)

for (b in 1:B) {
  
  cat("Iteration:", b, "\n")
  
  pop_b <- pop   # IMPORTANT: fresh copy each iteration
####################################
#### Stage 1: to Single disease ####
####################################
# From: OH node
n <- nrow(pop_b)

# Transition 1
# GenGamma: log(Ti)=μ0+Xiβ+σZi
mu0 <- 3.877369 
sigma <- 1.40459
Q <- 0.158643
beta_age <- -0.05655
beta_sex<- -0.09307
beta_wimd1<- -0.09795
beta_wimd2<- -0.02379
beta_wimd4<- 0.047789
beta_wimd5<- 0.218979

pop_b$lp1 <- mu0 + beta_age*pop_b$agec + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp1_exp<- exp(pop_b$lp1)
pop_b$t1 <- rgengamma(n = n, mu = pop_b$lp1, sigma = sigma, Q=Q)

# Transition 2
# Log-logistic: log(Ti)=μ0+Xiβ+1/γ*log(U/1−U)
scale <- 170.1529 #scale
mu0<- log(scale)
shape <- 0.877298 #shape
beta_age <- -0.04204
beta_sex<- 0.460691
beta_wimd1<- -0.33947
beta_wimd2<- -0.17021
beta_wimd4<- 0.103126
beta_wimd5<- 0.299392


pop_b$lp2 <- mu0 + beta_age*pop_b$agec + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp2_exp<-exp(pop_b$lp2)
pop_b$t2 <- flexsurv::rllogis(n = n, shape = shape, scale = pop_b$lp2_exp)

# Transition 3
# Log-logistic: log(Ti)=μ0+Xiβ+1/γ*log(U/1−U)
# GenGamma: log(Ti)=μ0+Xiβ+σZi
mu0 <- 5.312612 
sigma <- 1.749995
Q <- 0.513518
beta_age <- -0.08515
beta_sex<- 0.272441
beta_wimd1<- -0.07624
beta_wimd2<- -0.04557
beta_wimd4<- 0.017153
beta_wimd5<- 0.080824

pop_b$lp3 <- mu0 + beta_age*pop_b$agec + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp3_exp<- exp(pop_b$lp3)
pop_b$t3 <- rgengamma(n = n, mu = pop_b$lp3, sigma = sigma, Q=Q)

# Transition 4
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.096179 #shape
lambda <- 0.025202 #rate
beta_age <- 0.111307
beta_sex<- -0.36035
beta_wimd1<- 0.253223
beta_wimd2<- 0.106261
beta_wimd4<- -0.08651
beta_wimd5<- -0.20852

pop_b$lp4 <- log(lambda) + beta_age*pop_b$agec + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp4_exp <- exp(pop_b$lp4)
pop_b$t4<- rgompertz(n = n, shape = gamma, rate = pop_b$lp4_exp)

# Minimum path

tmat<-cbind(pop_b$t1,pop_b$t2,pop_b$t3,pop_b$t4)
tmat.min<-apply(tmat,1,min)
tmat.min.path<-apply(tmat,1,which.min)
tabulate(tmat.min.path)
pop_b$time1<-tmat.min
pop_b$trans1<- tmat.min.path
pop_b = subset(pop_b, select = -c(lp1,lp1_exp,lp2,lp2_exp,lp3,lp3_exp,lp4,lp4_exp) )

####################################
#### Stage 2: to Double disease ####
####################################
# Do not subset here? Simulate for everyone regardless of path taken and sort out path later?

#From: CKD node
pop_b$age2<- pop_b$age+pop_b$time1
pop_b$agec2<- pop_b$age2-65

# Transition 5
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- -0.03016 #shape
lambda <- 0.013595 #rate
beta_age <- 0.001276
beta_sex<- -0.34379
beta_wimd1<- 0.288066
beta_wimd2<- 0.169115
beta_wimd4<- -0.12399
beta_wimd5<- -0.11847


pop_b$lp5 <- log(lambda) + beta_age*pop_b$agec2 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp5_exp <- exp(pop_b$lp5)
pop_b$t5<- rgompertz(n = n, shape = gamma, rate = pop_b$lp5_exp)


# Transition 6
# GenGamma: log(Ti)=μ0+Xiβ+σZi
mu0 <- 4.657185 
sigma <- 1.526628
Q <- 0.55109
beta_age <- -0.05653
beta_sex<- 0.263712
beta_wimd1<- -0.02282
beta_wimd2<- 0.046647
beta_wimd4<- 0.10815
beta_wimd5<- 0.065599


pop_b$lp6 <- mu0 + beta_age*pop_b$agec2 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp6_exp<- exp(pop_b$lp6)
pop_b$t6 <- rgengamma(n = n, mu = pop_b$lp6, sigma = sigma, Q=Q)


# Transition 7
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.241107 #shape
lambda <- 0.014075 #rate
beta_age <- 0.125901
beta_sex<- -0.44279
beta_wimd1<- 0.179558
beta_wimd2<- 0.050926
beta_wimd4<- -0.04664
beta_wimd5<- -0.09087

pop_b$lp7 <- log(lambda) + beta_age*pop_b$agec2 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp7_exp <- exp(pop_b$lp7)
pop_b$t7<- rgompertz(n = n, shape = gamma, rate = pop_b$lp7_exp)

pop_b = subset(pop_b, select = -c(lp5,lp5_exp,lp6,lp6_exp,lp7,lp7_exp) )

#From: D node

# Transition 8
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- -0.13642 #shape
lambda <- 0.037056 #rate
beta_age <- 0.049365
beta_sex<- 0.187662
beta_wimd1<- 0.079974
beta_wimd2<- 0.053256
beta_wimd4<- -0.0009
beta_wimd5<- -0.14991

pop_b$lp8 <- log(lambda) + beta_age*pop_b$agec2 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp8_exp <- exp(pop_b$lp8)
pop_b$t8<- rgompertz(n = n, shape = gamma, rate = pop_b$lp8_exp)

# Transition 9
# Log-logistic: log(Ti)=μ0+Xiβ+1/γ*log(U/1−U)
scale <- 68.37181 #scale
mu0<- log(scale)
shape <- 0.997204 #shape
beta_age <- -0.05896
beta_sex<- 0.188338
beta_wimd1<- -0.04488
beta_wimd2<- -0.05505
beta_wimd4<- -0.004
beta_wimd5<- 0.140477


pop_b$lp9 <- mu0 + beta_age*pop_b$agec + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp9_exp<-exp(pop_b$lp9)
pop_b$t9 <- flexsurv::rllogis(n = n, shape = shape, scale = pop_b$lp9_exp)


# Transition 10
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.140436 #shape
lambda <- 0.026176 #rate
beta_age <- 0.115412
beta_sex<- -0.2688
beta_wimd1<- 0.254
beta_wimd2<- 0.037903
beta_wimd4<- -0.15311
beta_wimd5<- -0.23841

pop_b$lp10 <- log(lambda) + beta_age*pop_b$agec2 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp10_exp <- exp(pop_b$lp10)
pop_b$t10<- rgompertz(n = n, shape = gamma, rate = pop_b$lp10_exp)

pop_b = subset(pop_b, select = -c(lp8,lp8_exp,lp9,lp9_exp,lp10,lp10_exp) )


#From: HF node

# Transition 11
# Log-logistic: log(Ti)=μ0+Xiβ+1/γ*log(U/1−U)
scale <- 21.19925 #scale
mu0<- log(scale)
shape <- 1.166593 #shape
beta_age <- -0.0468
beta_sex<- 0.005295
beta_wimd1<- 0.004856
beta_wimd2<- 0.068333
beta_wimd4<- 0.056276
beta_wimd5<- 0.133569

pop_b$lp11 <- mu0 + beta_age*pop_b$agec + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp11_exp<-exp(pop_b$lp11)
pop_b$t11 <- flexsurv::rllogis(n = n, shape = shape, scale = pop_b$lp11_exp)

# Transition 12
# Log-logistic: log(Ti)=μ0+Xiβ+1/γ*log(U/1−U)
scale <- 88.88381 #scale
mu0<- log(scale)
shape <- 0.970517 #shape
beta_age <- -0.00329
beta_sex<- 0.258128
beta_wimd1<- -0.24214
beta_wimd2<- -0.16556
beta_wimd4<- 0.054737
beta_wimd5<- 0.103101

pop_b$lp12 <- mu0 + beta_age*pop_b$agec + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp12_exp<-exp(pop_b$lp12)
pop_b$t12 <- flexsurv::rllogis(n = n, shape = shape, scale = pop_b$lp12_exp)

# Transition 13
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.202875 #shape
lambda <- 0.037837 #rate
beta_age <- 0.092137
beta_sex<- -0.26984
beta_wimd1<- 0.176944
beta_wimd2<- 0.050968
beta_wimd4<- -0.08867
beta_wimd5<- -0.19682

pop_b$lp13 <- log(lambda) + beta_age*pop_b$agec2 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp13_exp <- exp(pop_b$lp13)
pop_b$t13<- rgompertz(n = n, shape = gamma, rate = pop_b$lp13_exp)

pop_b = subset(pop_b, select = -c(lp11,lp11_exp,lp12,lp12_exp,lp13,lp13_exp))

# matrices
m1 <- as.matrix(pop_b[, c("t5","t6","t7")])
m2 <- as.matrix(pop_b[, c("t8","t9","t10")])
m3 <- as.matrix(pop_b[, c("t11","t12","t13")])

# row-wise minima
min1 <- do.call(pmin, c(as.data.frame(m1), na.rm = TRUE))
min2 <- do.call(pmin, c(as.data.frame(m2), na.rm = TRUE))
min3 <- do.call(pmin, c(as.data.frame(m3), na.rm = TRUE))

# indices
idx1 <- max.col(-m1, ties.method = "first")  # 1–3
idx2 <- max.col(-m2, ties.method = "first")  # 1–3
idx3 <- max.col(-m3, ties.method = "first")  # 1–3

# handle all-NA rows
all_na1 <- rowSums(!is.na(m1)) == 0
all_na2 <- rowSums(!is.na(m2)) == 0
all_na3 <- rowSums(!is.na(m3)) == 0

min1[all_na1] <- NA
min2[all_na2] <- NA
min3[all_na3] <- NA

idx1[all_na1] <- NA
idx2[all_na2] <- NA
idx3[all_na3] <- NA

pop_b$time2 <- ifelse(pop_b$trans1 == 1, min1,
                    ifelse(pop_b$trans1 == 2, min2,
                           ifelse(pop_b$trans1 == 3, min3, NA)))

pop_b$trans2 <- ifelse(pop_b$trans1 == 1, c(5L,6L,7L)[idx1],
                     ifelse(pop_b$trans1 == 2, c(8L,9L,10L)[idx2],
                            ifelse(pop_b$trans1 == 3, c(11L,12L,13L)[idx3], NA)))

pop_b$time2_age <- pop_b$time2
pop_b$time2_age[is.na(pop_b$time2)]<- 0

####################################
#### Stage 3: to Triple disease ####
####################################

pop_b$age3 <- pop_b$age2 + pop_b$time2_age
pop_b$agec3 <- pop_b$age3-65

# From: CKD + D node

# Transition 14
# Exponential: h(t∣X)=h0(t)exp(Xβ)
rate<- 0.014601
beta_age <- 0.049312
beta_sex<- 0.003611
beta_wimd1<- -0.25765
beta_wimd2<- -0.52523
beta_wimd4<- -0.38659
beta_wimd5<- -0.08755

pop_b$lp14<- rate + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp14_exp <- exp(pop_b$lp14)

pop_b$t14<-rexp(n,rate = pop_b$lp14_exp)

# Transition 15
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.195135 #shape
lambda <- 0.02186 #rate
beta_age <- 0.125166
beta_sex<- -0.29602
beta_wimd1<- 0.223866
beta_wimd2<- -0.05843
beta_wimd4<- -0.08135
beta_wimd5<- 0.035473

pop_b$lp15 <- log(lambda) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp15_exp <- exp(pop_b$lp15)
pop_b$t15<- rgompertz(n = n, shape = gamma, rate = pop_b$lp15_exp)


# From: CKD + HF node

# Transition 16
# Exponential: h(t∣X)=h0(t)exp(Xβ)
rate<- 0.018451
beta_age <- -0.01669
beta_sex<- -0.3062
beta_wimd1<- 0.03038
beta_wimd2<- 0.213786
beta_wimd4<- -0.35421
beta_wimd5<- -0.03037

pop_b$lp16<- rate + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp16_exp <- exp(pop_b$lp16)
pop_b$t16<-rexp(n,rate = pop_b$lp16_exp)

# Transition 17
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.242813 #shape
lambda <- 0.057286 #rate
beta_age <- 0.089584
beta_sex<- -0.39712
beta_wimd1<- 0.067646
beta_wimd2<- -0.02114
beta_wimd4<- -0.19859
beta_wimd5<- -0.2892

pop_b$lp17 <- log(lambda) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp17_exp <- exp(pop_b$lp17)
pop_b$t17<- rgompertz(n = n, shape = gamma, rate = pop_b$lp17_exp)


# From: D + CKD node

# Transition 18
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- -0.07434 #shape
lambda <- 0.014296 #rate
beta_age <- 0.048322
beta_sex<- -0.07746
beta_wimd1<- 0.212806
beta_wimd2<- 0.32492
beta_wimd4<- 0.379758
beta_wimd5<- 0.183934

pop_b$lp18 <- log(lambda) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp18_exp <- exp(pop_b$lp18)
pop_b$t18<- rgompertz(n = n, shape = gamma, rate = pop_b$lp18_exp)

# Transition 19
# GenGamma: log(Ti)=μ0+Xiβ+σZi
mu0 <- 2.725085 
sigma <- 0.577348
Q <- 1.282546
beta_age <- -0.07936
beta_sex<- 0.206411
beta_wimd1<- 0.009628
beta_wimd2<- -0.00885
beta_wimd4<- 0.104093
beta_wimd5<- 0.163409

pop_b$lp19 <- mu0 + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp19_exp<- exp(pop_b$lp19)
pop_b$t19 <- rgengamma(n = n, mu = pop_b$lp19, sigma = sigma, Q=Q)


# From: D + HF node

# Transition 20
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- -0.16187 #shape
lambda <- 0.047619 #rate
beta_age <- 0.057682
beta_sex<- 0.10327
beta_wimd1<- -0.00077
beta_wimd2<- 0.275854
beta_wimd4<- 0.115321
beta_wimd5<- -0.53076

pop_b$lp20 <- log(lambda) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp20_exp <- exp(pop_b$lp20)
pop_b$t20<- rgompertz(n = n, shape = gamma, rate = pop_b$lp20_exp)

# Transition 21
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.205247 #shape
lambda <- 0.051003 #rate
beta_age <- 0.092699
beta_sex<- -0.27207
beta_wimd1<- 0.144285
beta_wimd2<- 0.171476
beta_wimd4<- -0.29013
beta_wimd5<- -0.26129

pop_b$lp21 <- log(lambda) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp21_exp <- exp(pop_b$lp21)
pop_b$t21<- rgompertz(n = n, shape = gamma, rate = pop_b$lp21_exp)



# From: HF + CKD node

# Transition 22
# Exponential: h(t∣X)=h0(t)exp(Xβ)
rate<- 0.015979
beta_age <- -0.00804
beta_sex<- -0.15176
beta_wimd1<- 0.484209
beta_wimd2<- 0.174044
beta_wimd4<- -0.13086
beta_wimd5<- 0.051782

pop_b$lp22<- rate + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp22_exp <- exp(pop_b$lp22)
pop_b$t22<-rexp(n,rate = pop_b$lp22_exp)

# Transition 23
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.321955 #shape
lambda <- 0.038839 #rate
beta_age <- 0.076558
beta_sex<- -0.23807
beta_wimd1<- 0.152193
beta_wimd2<- 0.041942
beta_wimd4<- -0.09044
beta_wimd5<- -0.08103

pop_b$lp23 <- log(lambda) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp23_exp <- exp(pop_b$lp23)
pop_b$t23<- rgompertz(n = n, shape = gamma, rate = pop_b$lp23_exp)


# From: HF + D node

# Transition 24
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- -0.18247 #shape
lambda <- 0.060587 #rate
beta_age <- 0.04846
beta_sex<- -0.09198
beta_wimd1<- -0.17006
beta_wimd2<- -0.11631
beta_wimd4<- -0.0839
beta_wimd5<- -0.17923

pop_b$lp24 <- log(lambda) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp24_exp <- exp(pop_b$lp24)
pop_b$t24<- rgompertz(n = n, shape = gamma, rate = pop_b$lp24_exp)

# Transition 25
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.176138 #shape
lambda <- 0.067468 #rate
beta_age <- 0.078933
beta_sex<- -0.18335
beta_wimd1<- 0.112176
beta_wimd2<- 0.154216
beta_wimd4<- 0.04644
beta_wimd5<- -0.19369

pop_b$lp25 <- log(lambda) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp25_exp <- exp(pop_b$lp25)
pop_b$t25<- rgompertz(n = n, shape = gamma, rate = pop_b$lp25_exp)

pop_b = subset(pop_b, select = -c(lp14,lp14_exp,lp15,lp15_exp,lp16,lp16_exp,lp17,lp17_exp,lp18,lp18_exp,lp19,lp19_exp,lp20,lp20_exp,lp21,lp21_exp,lp22,lp22_exp,lp23,lp23_exp,lp24,lp24_exp,lp25,lp25_exp))

# matrices
m1 <- as.matrix(pop_b[, c("t14","t15")])
m2 <- as.matrix(pop_b[, c("t16","t17")])
m3 <- as.matrix(pop_b[, c("t18","t19")])
m4 <- as.matrix(pop_b[, c("t20","t21")])
m5 <- as.matrix(pop_b[, c("t22","t23")])
m6 <- as.matrix(pop_b[, c("t24","t25")])


# row-wise minima
min1 <- do.call(pmin, c(as.data.frame(m1), na.rm = TRUE))
min2 <- do.call(pmin, c(as.data.frame(m2), na.rm = TRUE))
min3 <- do.call(pmin, c(as.data.frame(m3), na.rm = TRUE))
min4 <- do.call(pmin, c(as.data.frame(m4), na.rm = TRUE))
min5 <- do.call(pmin, c(as.data.frame(m5), na.rm = TRUE))
min6 <- do.call(pmin, c(as.data.frame(m6), na.rm = TRUE))

# indices
idx1 <- max.col(-m1, ties.method = "first")  # 1–3
idx2 <- max.col(-m2, ties.method = "first")  # 1–3
idx3 <- max.col(-m3, ties.method = "first")  # 1–3
idx4 <- max.col(-m4, ties.method = "first")  # 1–3
idx5 <- max.col(-m5, ties.method = "first")  # 1–3
idx6 <- max.col(-m6, ties.method = "first")  # 1–3

# handle all-NA rows
all_na1 <- rowSums(!is.na(m1)) == 0
all_na2 <- rowSums(!is.na(m2)) == 0
all_na3 <- rowSums(!is.na(m3)) == 0
all_na4 <- rowSums(!is.na(m4)) == 0
all_na5 <- rowSums(!is.na(m5)) == 0
all_na6 <- rowSums(!is.na(m6)) == 0

min1[all_na1] <- NA
min2[all_na2] <- NA
min3[all_na3] <- NA
min4[all_na4] <- NA
min5[all_na5] <- NA
min6[all_na6] <- NA

idx1[all_na1] <- NA
idx2[all_na2] <- NA
idx3[all_na3] <- NA
idx4[all_na4] <- NA
idx5[all_na5] <- NA
idx6[all_na6] <- NA

pop_b$time3 <- ifelse(pop_b$trans2 == 5, min1,
                    ifelse(pop_b$trans2 == 6, min2, 
                           ifelse(pop_b$trans2==8, min3,
                                  ifelse(pop_b$trans2==9, min4,
                                         ifelse(pop_b$trans2==11, min5,
                                                ifelse(pop_b$trans2==12,min6,
                           NA))))))

pop_b$trans3 <- ifelse(pop_b$trans2 == 5, c(14L,15L)[idx1],
                     ifelse(pop_b$trans2 == 6, c(16L,17L)[idx2],
                            ifelse(pop_b$trans2==8, c(18L,19L)[idx3],
                                   ifelse(pop_b$trans2==9, c(20L,21L)[idx4],
                                          ifelse(pop_b$trans2==11, c(22L,23L)[idx5],
                                                 ifelse(pop_b$trans2==12, c(24L,25L)[idx6],
                                          NA))))))


pop_b$time3_age <- pop_b$time3
pop_b$time3_age[is.na(pop_b$time3)]<- 0

###########################
#### Stage 4: to death ####
###########################

pop_b$age4 <- pop_b$age3 + pop_b$time3_age
pop_b$agec4 <- pop_b$age4-65


# From: CKD + D + HF node
# Transition 26
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.268892 #shape
lambda <- 0.059648 #rate
beta_age <- 0.079816
beta_sex<- 0.059648
beta_wimd1<- 0.638133
beta_wimd2<- -0.31525
beta_wimd4<- 0.119112
beta_wimd5<- 0.434409

pop_b$lp26 <- log(lambda) + beta_age*pop_b$agec4 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp26_exp <- exp(pop_b$lp26)
pop_b$t26<- rgompertz(n = n, shape = gamma, rate = pop_b$lp26_exp)


# From: CKD + HF + D node
# Transition 27
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.280503 #shape
lambda <- 0.084889 #rate
beta_age <- 0.069536
beta_sex<- -0.00551
beta_wimd1<- 0.4453
beta_wimd2<- 0.119993
beta_wimd4<- -0.25663
beta_wimd5<- -0.3128

pop_b$lp27 <- log(lambda) + beta_age*pop_b$agec4 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp27_exp <- exp(pop_b$lp27)
pop_b$t27<- rgompertz(n = n, shape = gamma, rate = pop_b$lp27_exp)


# From: D + CKD + HF node
# Transition 28
# GenGamma: log(Ti)=μ0+Xiβ+σZi
mu0 <- 1.96282 
sigma <- 0.528771
Q <- 1.938412
beta_age <- -0.04778
beta_sex<- 0.015592
beta_wimd1<- 0.187682
beta_wimd2<- 0.142381
beta_wimd4<- 0.131647
beta_wimd5<- -0.13457

pop_b$lp28 <- mu0 + beta_age*pop_b$agec4 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp28_exp<- exp(pop_b$lp28)
pop_b$t28 <- rgengamma(n = n, mu = pop_b$lp28, sigma = sigma, Q=Q)

# From: D + HF + CKD node

# Transition 29
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.429129 #shape
lambda <- 0.045404 #rate
beta_age <- 0.074555
beta_sex<- -0.44637
beta_wimd1<- 0.067292
beta_wimd2<- -0.00642
beta_wimd4<- -0.23454
beta_wimd5<- -0.54576

pop_b$lp29 <- log(lambda) + beta_age*pop_b$agec4 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp29_exp <- exp(pop_b$lp29)
pop_b$t29<- rgompertz(n = n, shape = gamma, rate = pop_b$lp29_exp)


# From: HF + CKD + D node
# Transition 30
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.349587 #shape
lambda <- 0.069599 #rate
beta_age <- 0.092623
beta_sex<- -0.42267
beta_wimd1<- -0.10068
beta_wimd2<- -0.00528
beta_wimd4<- 0.128351
beta_wimd5<- 0.035094

pop_b$lp30 <- log(lambda) + beta_age*pop_b$agec4 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp30_exp <- exp(pop_b$lp30)
pop_b$t30<- rgompertz(n = n, shape = gamma, rate = pop_b$lp30_exp)


# From HF + D + CKD
# Transition 31
# Gompertz: t=1/γ*log(1−γ/λ*(log(U)/exp(Xβ))
gamma<- 0.357486 #shape
lambda <- 0.042954 #rate
beta_age <- 0.074485
beta_sex<- -0.37109
beta_wimd1<- 0.299921
beta_wimd2<- -0.13012
beta_wimd4<- 0.253627
beta_wimd5<- 0.424614

pop_b$lp31 <- log(lambda) + beta_age*pop_b$agec4 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
pop_b$lp31_exp <- exp(pop_b$lp31)
pop_b$t31<- rgompertz(n = n, shape = gamma, rate = pop_b$lp31_exp)

pop_b = subset(pop_b, select = -c(lp26,lp26_exp,lp27,lp27_exp,lp28,lp28_exp,lp29,lp29_exp,lp30,lp30_exp,lp31,lp31_exp))

pop_b <- pop_b %>%
  mutate(
    time4 = case_when(
      trans3 == 14 ~ t26,
      trans3 == 16 ~ t27,
      trans3 == 18 ~ t28,
      trans3 == 20 ~ t29,
      trans3 == 22 ~ t30,
      trans3 == 24 ~ t31,
      TRUE ~ NA_real_
    )
  )

pop_b <- pop_b %>%
  mutate(
    trans4 = case_when(
      trans3 == 14 ~ 26,
      trans3 == 16 ~ 27,
      trans3 == 18 ~ 28,
      trans3 == 20 ~ 29,
      trans3 == 22 ~ 30,
      trans3 == 24 ~ 31,
      TRUE ~ NA_real_
    )
  )

pop_b = subset(pop_b, select = -c(agec,agec2,agec3,agec4,time2_age,time3_age))

pop_b$totly<- rowSums(pop_b[,c("age","time1", "time2", "time3", "time4")], na.rm=TRUE)
pop_b$life.exp<-pop_b$totly

# Results for untreated
trans_levels2 <- c(5,6,7,8,9,10,11,12,13)
tab2 <- table(factor(pop_b$trans2, levels = trans_levels2))
prop2 <- as.numeric(prop.table(tab2))

trans_levels3 <- c(14,15,16,17,18,19,20,21,22,23,24,25)
tab3 <- table(factor(pop_b$trans3, levels = trans_levels3))
prop3 <- as.numeric(prop.table(tab3))

trans_levels4 <- c(26,27,28,29,30,31)
tab4 <- table(factor(pop_b$trans4, levels = trans_levels4))
prop4 <- as.numeric(prop.table(tab4))

results[[b]] <- data.frame(
  mean_life_exp = mean(pop_b$life.exp),
  # restricted to 100 years
  rmst = mean(pmin(pop_b$life.exp, 100)),
  mean_t5 = mean(pop_b$time2[pop_b$trans2==5], na.rm=TRUE),
  mean_t6 = mean(pop_b$time2[pop_b$trans2==6], na.rm=TRUE),
  mean_t7 = mean(pop_b$time2[pop_b$trans2==7], na.rm=TRUE),
  mean_t8 = mean(pop_b$time2[pop_b$trans2==8], na.rm=TRUE),
  mean_t9 = mean(pop_b$time2[pop_b$trans2==9], na.rm=TRUE),
  mean_t10 = mean(pop_b$time2[pop_b$trans2==10], na.rm=TRUE),
  mean_t11 = mean(pop_b$time2[pop_b$trans2==11], na.rm=TRUE),
  mean_t12 = mean(pop_b$time2[pop_b$trans2==12], na.rm=TRUE),
  mean_t13 = mean(pop_b$time2[pop_b$trans2==13], na.rm=TRUE),
  mean_t14 = mean(pop_b$time3[pop_b$trans3==14], na.rm=TRUE),
  mean_t15 = mean(pop_b$time3[pop_b$trans3==15], na.rm=TRUE),
  mean_t16 = mean(pop_b$time3[pop_b$trans3==16], na.rm=TRUE),
  mean_t17 = mean(pop_b$time3[pop_b$trans3==17], na.rm=TRUE),
  mean_t18 = mean(pop_b$time3[pop_b$trans3==18], na.rm=TRUE),
  mean_t19 = mean(pop_b$time3[pop_b$trans3==19], na.rm=TRUE),
  mean_t20 = mean(pop_b$time3[pop_b$trans3==20], na.rm=TRUE),
  mean_t21 = mean(pop_b$time3[pop_b$trans3==21], na.rm=TRUE),
  mean_t22 = mean(pop_b$time3[pop_b$trans3==22], na.rm=TRUE),
  mean_t23 = mean(pop_b$time3[pop_b$trans3==23], na.rm=TRUE),
  mean_t24 = mean(pop_b$time3[pop_b$trans3==24], na.rm=TRUE),
  mean_t25 = mean(pop_b$time3[pop_b$trans3==25], na.rm=TRUE),
  mean_t26 = mean(pop_b$time4[pop_b$trans4==26], na.rm=TRUE),
  mean_t27 = mean(pop_b$time4[pop_b$trans4==27], na.rm=TRUE),
  mean_t28 = mean(pop_b$time4[pop_b$trans4==28], na.rm=TRUE),
  mean_t29 = mean(pop_b$time4[pop_b$trans4==29], na.rm=TRUE),
  mean_t30 = mean(pop_b$time4[pop_b$trans4==30], na.rm=TRUE),
  mean_t31 = mean(pop_b$time4[pop_b$trans4==31], na.rm=TRUE),
  t5 = prop2[1],
  t6 = prop2[2],
  t7 = prop2[3],
  t8 = prop2[4],
  t9 = prop2[5],
  t10 = prop2[6],
  t11 = prop2[7],
  t12 = prop2[8],
  t13 = prop2[9],
  t14 = prop3[1],
  t15 = prop3[2],
  t16 = prop3[3],
  t17 = prop3[4],
  t18 = prop3[5], 
  t19 = prop3[6],
  t20 = prop3[7],
  t21 = prop3[8], 
  t22 = prop3[9],
  t23 = prop3[10],
  t24 = prop3[11],
  t25 = prop3[12], 
  t26 = prop4[1],
  t27 = prop4[2],
  t28 = prop4[3],
  t29 = prop4[4],
  t30 = prop4[5], 
  t31 = prop4[6]
)

mc_pop[[b]] <- pop_b[, c("time2", "trans2", "time3", "trans3", "time4", "trans4")]
}

res <- do.call(rbind, results)

summary(res)
apply(res, 2, quantile, c(0.025, 0.5, 0.975))


###################################

#####################################
########### Treated arm #############
#####################################
pop<-pop_b

hr_table <- data.frame(
  trans = 5:31,
  hr = c(
    0.67,0.71,0.69,0.76,0.83,0.93,0.71,
    0.67,0.83,0.78,0.91,0.67,0.56,
    0.78,0.91,0.58,0.87,0.67,0.56,
    0.58,0.87,0.91,0.91,0.91,0.91,
    0.91,0.91
  ),
  se_loghr = c(
    0.14,0.13,0.13,0.07,0.07,0.06,0.25,
    0.14,0.08,0.17,0.16,0.14,0.26,
    0.17,0.16,0.24,0.13,0.14,0.26,
    0.24,0.13,0.16,0.16,0.16,0.16,
    0.16,0.16
  )
)

# Monte Carlo wrapper #
B <- 100   # number of simulations
n <- nrow(pop)

results_trt <- vector("list", B)
mc_pop_trt  <- vector("list", B)

sample_HR <- function(hr_table) {
  
  log_hr <- log(hr_table$hr)
  sampled_log_hr <- rnorm(nrow(hr_table), mean = log_hr, sd = hr_table$se_loghr)
  
  exp(sampled_log_hr)
}

for (b in 1:B) {
  
  cat("Iteration:", b, "\n")
  
  pop_b <- pop   # IMPORTANT: fresh copy each iteration
  
  HR_b <- sample_HR(hr_table)
  
  # --------------------------
  # STAGE 1–4
  # --------------------------
  # BUT replace each HR with HR_b[trans] to allow for uncertainty
  # --------------------------
  
  ## Example replacement pattern:
  ## instead of HR <- 0.67
  ## use HR <- HR_b[5]  (for transition 5)

  HR <- HR_b[1]
  sf <- survfit(Surv(t5, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t5_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[2]
  sf <- survfit(Surv(t6, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t6_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[3]
  sf <- survfit(Surv(t7, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t7_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[4]
  sf <- survfit(Surv(t8, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t8_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[5]
  sf <- survfit(Surv(t9, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t9_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[6]
  sf <- survfit(Surv(t10, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t10_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[7]
  sf <- survfit(Surv(t11, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t11_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[8]
  sf <- survfit(Surv(t12, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t12_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[9]
  sf <- survfit(Surv(t13, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t13_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[10]
  sf <- survfit(Surv(t14, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t14_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[11]
  sf <- survfit(Surv(t15, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t15_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[12]
  sf <- survfit(Surv(t16, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t16_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[13]
  sf <- survfit(Surv(t17, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t17_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[14]
  sf <- survfit(Surv(t18, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t18_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[15]
  sf <- survfit(Surv(t19, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t19_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[16]
  sf <- survfit(Surv(t20, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t20_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[17]
  sf <- survfit(Surv(t21, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t21_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[18]
  sf <- survfit(Surv(t22, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t22_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[19]
  sf <- survfit(Surv(t23, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t23_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[20]
  sf <- survfit(Surv(t24, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t24_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))

  HR <- HR_b[21]
  sf <- survfit(Surv(t25, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t25_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[22]
  sf <- survfit(Surv(t26, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t26_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[23]
  sf <- survfit(Surv(t27, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t27_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[24]
  sf <- survfit(Surv(t28, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t28_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[25]
  sf <- survfit(Surv(t29, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t29_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[26]
  sf <- survfit(Surv(t30, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t30_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  HR <- HR_b[27]
  sf <- survfit(Surv(t31, rep(1,n)) ~ 1, data = pop_b)
  pop_b$t31_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(runif(n))
  
  # matrices
  m1 <- as.matrix(pop_b[, c("t5_trt","t6_trt","t7_trt")])
  m2 <- as.matrix(pop_b[, c("t8_trt","t9_trt","t10_trt")])
  m3 <- as.matrix(pop_b[, c("t11_trt","t12_trt","t13_trt")])
  
  # row-wise minima
  min1 <- do.call(pmin, c(as.data.frame(m1), na.rm = TRUE))
  min2 <- do.call(pmin, c(as.data.frame(m2), na.rm = TRUE))
  min3 <- do.call(pmin, c(as.data.frame(m3), na.rm = TRUE))
  
  # indices
  idx1 <- max.col(-m1, ties.method = "first")  # 1–3
  idx2 <- max.col(-m2, ties.method = "first")  # 1–3
  idx3 <- max.col(-m3, ties.method = "first")  # 1–3
  
  # handle all-NA rows
  all_na1 <- rowSums(!is.na(m1)) == 0
  all_na2 <- rowSums(!is.na(m2)) == 0
  all_na3 <- rowSums(!is.na(m3)) == 0
  
  min1[all_na1] <- NA
  min2[all_na2] <- NA
  min3[all_na3] <- NA
  
  idx1[all_na1] <- NA
  idx2[all_na2] <- NA
  idx3[all_na3] <- NA
  
  pop_b$time2_trt <- ifelse(pop_b$trans1 == 1, min1,
                            ifelse(pop_b$trans1 == 2, min2,
                                   ifelse(pop_b$trans1 == 3, min3, NA)))
  
  pop_b$trans2_trt <- ifelse(pop_b$trans1 == 1, c(5L,6L,7L)[idx1],
                             ifelse(pop_b$trans1 == 2, c(8L,9L,10L)[idx2],
                                    ifelse(pop_b$trans1 == 3, c(11L,12L,13L)[idx3], NA)))
  
  # matrices
  m1 <- as.matrix(pop_b[, c("t14_trt","t15_trt")])
  m2 <- as.matrix(pop_b[, c("t16_trt","t17_trt")])
  m3 <- as.matrix(pop_b[, c("t18_trt","t19_trt")])
  m4 <- as.matrix(pop_b[, c("t20_trt","t21_trt")])
  m5 <- as.matrix(pop_b[, c("t22_trt","t23_trt")])
  m6 <- as.matrix(pop_b[, c("t24_trt","t25_trt")])
  
  
  # row-wise minima
  min1 <- do.call(pmin, c(as.data.frame(m1), na.rm = TRUE))
  min2 <- do.call(pmin, c(as.data.frame(m2), na.rm = TRUE))
  min3 <- do.call(pmin, c(as.data.frame(m3), na.rm = TRUE))
  min4 <- do.call(pmin, c(as.data.frame(m4), na.rm = TRUE))
  min5 <- do.call(pmin, c(as.data.frame(m5), na.rm = TRUE))
  min6 <- do.call(pmin, c(as.data.frame(m6), na.rm = TRUE))
  
  # indices
  idx1 <- max.col(-m1, ties.method = "first")  # 1–3
  idx2 <- max.col(-m2, ties.method = "first")  # 1–3
  idx3 <- max.col(-m3, ties.method = "first")  # 1–3
  idx4 <- max.col(-m4, ties.method = "first")  # 1–3
  idx5 <- max.col(-m5, ties.method = "first")  # 1–3
  idx6 <- max.col(-m6, ties.method = "first")  # 1–3
  
  # handle all-NA rows
  all_na1 <- rowSums(!is.na(m1)) == 0
  all_na2 <- rowSums(!is.na(m2)) == 0
  all_na3 <- rowSums(!is.na(m3)) == 0
  all_na4 <- rowSums(!is.na(m4)) == 0
  all_na5 <- rowSums(!is.na(m5)) == 0
  all_na6 <- rowSums(!is.na(m6)) == 0
  
  min1[all_na1] <- NA
  min2[all_na2] <- NA
  min3[all_na3] <- NA
  min4[all_na4] <- NA
  min5[all_na5] <- NA
  min6[all_na6] <- NA
  
  idx1[all_na1] <- NA
  idx2[all_na2] <- NA
  idx3[all_na3] <- NA
  idx4[all_na4] <- NA
  idx5[all_na5] <- NA
  idx6[all_na6] <- NA
  
  pop_b$time3_trt <- ifelse(pop_b$trans2_trt == 5, min1,
                          ifelse(pop_b$trans2_trt == 6, min2, 
                                 ifelse(pop_b$trans2_trt ==8, min3,
                                        ifelse(pop_b$trans2_trt ==9, min4,
                                               ifelse(pop_b$trans2_trt ==11, min5,
                                                      ifelse(pop_b$trans2_trt ==12,min6,
                                                             NA))))))
  
  pop_b$trans3_trt <- ifelse(pop_b$trans2_trt == 5, c(14L,15L)[idx1],
                           ifelse(pop_b$trans2_trt == 6, c(16L,17L)[idx2],
                                  ifelse(pop_b$trans2_trt ==8, c(18L,19L)[idx3],
                                         ifelse(pop_b$trans2_trt ==9, c(20L,21L)[idx4],
                                                ifelse(pop_b$trans2_trt ==11, c(22L,23L)[idx5],
                                                       ifelse(pop_b$trans2_trt ==12, c(24L,25L)[idx6],
                                                              NA))))))
  
  
  pop_b <- pop_b %>%
    mutate(
      time4_trt = case_when(
        trans3_trt == 14 ~ t26_trt,
        trans3_trt == 16 ~ t27_trt,
        trans3_trt == 18 ~ t28_trt,
        trans3_trt == 20 ~ t29_trt,
        trans3_trt == 22 ~ t30_trt,
        trans3_trt == 24 ~ t31_trt,
        TRUE ~ NA_real_
      )
    )
  
  pop_b <- pop_b %>%
    mutate(
      trans4_trt = case_when(
        trans3_trt == 14 ~ 26,
        trans3_trt == 16 ~ 27,
        trans3_trt == 18 ~ 28,
        trans3_trt == 20 ~ 29,
        trans3_trt == 22 ~ 30,
        trans3_trt == 24 ~ 31,
        TRUE ~ NA_real_
      )
    )
  
  
  pop_b$totly_trt<- rowSums(pop_b[,c("age","time1", "time2_trt", "time3_trt", "time4_trt")], na.rm=TRUE)
  pop_b$life.exp_trt<-pop_b$totly_trt
  
  # Function to format estimate and CI
  fmt <- function(est, lci, uci, digits = 2) {
    sprintf("%.*f (%.*f, %.*f)", digits, est, digits, lci, digits, uci)
  }
  
  # Reshape untreated data
  untreated_summary <- pop_b %>%
    select(trans2, trans3, trans4, life.exp, age2, age3, age4) %>%
    pivot_longer(
      cols = c(trans2, trans3, trans4),
      names_to = "transition",
      values_to = "trans_value"
    ) %>%
    pivot_longer(
      cols = c(age2, age3, age4),
      names_to = "age_transition",
      values_to = "age"
    ) %>%
    filter(gsub("age", "trans", age_transition) == transition) %>%
    filter(
      (transition == "trans2" & trans_value >= 5  & trans_value <= 13) |
        (transition == "trans3" & trans_value >= 14 & trans_value <= 25) |
        (transition == "trans4" & trans_value >= 26 & trans_value <= 31)
    ) %>%
    group_by(transition, trans_value) %>%
    summarise(
      mean_life_exp = mean(life.exp, na.rm = TRUE),
      ci_lower_life_exp = mean(life.exp, na.rm = TRUE) - qt(0.975, df = n() - 1) * sd(life.exp, na.rm = TRUE) / sqrt(n()),
      ci_upper_life_exp = mean(life.exp, na.rm = TRUE) + qt(0.975, df = n() - 1) * sd(life.exp, na.rm = TRUE) / sqrt(n()),
      mean_age = mean(age, na.rm = TRUE),
      ci_lower_age = mean(age, na.rm = TRUE) - qt(0.975, df = n() - 1) * sd(age, na.rm = TRUE) / sqrt(n()),
      ci_upper_age = mean(age, na.rm = TRUE) + qt(0.975, df = n() - 1) * sd(age, na.rm = TRUE) / sqrt(n()),
      n = n(), .groups = "drop") %>%
    group_by(transition) %>%
    mutate(
      proportion = n / sum(n),
      ci_lower_prop = (2 * n + 1.96^2 - 1.96 * sqrt(1.96^2 + 4 * n * (1 - proportion))) / (2 * (sum(n) + 1.96^2)),
      ci_upper_prop = (2 * n + 1.96^2 + 1.96 * sqrt(1.96^2 + 4 * n * (1 - proportion))) / (2 * (sum(n) + 1.96^2))
    ) %>%
    ungroup() %>%
    mutate(group = "Untreated")
  
  # Reshape treated data
  treated_summary <- pop_b %>%
    select(trans2_trt, trans3_trt, trans4_trt, life.exp_trt, age2, age3, age4) %>%
    pivot_longer(
      cols = c(trans2_trt, trans3_trt, trans4_trt),
      names_to = "transition",
      values_to = "trans_value"
    ) %>%
    mutate(transition = gsub("_trt", "", transition)) %>%
    pivot_longer(
      cols = c(age2, age3, age4),
      names_to = "age_transition",
      values_to = "age"
    ) %>%
    filter(gsub("age", "trans", age_transition) == transition) %>%
    filter(
      (transition == "trans2" & trans_value >= 5  & trans_value <= 13) |
        (transition == "trans3" & trans_value >= 14 & trans_value <= 25) |
        (transition == "trans4" & trans_value >= 26 & trans_value <= 31)
    ) %>%
    group_by(transition, trans_value) %>%
    summarise(
      mean_life_exp = mean(life.exp_trt, na.rm = TRUE),
      ci_lower_life_exp = mean(life.exp_trt, na.rm = TRUE) - qt(0.975, df = n() - 1) * sd(life.exp_trt, na.rm = TRUE) / sqrt(n()),
      ci_upper_life_exp = mean(life.exp_trt, na.rm = TRUE) + qt(0.975, df = n() - 1) * sd(life.exp_trt, na.rm = TRUE) / sqrt(n()),
      mean_age = mean(age, na.rm = TRUE),
      ci_lower_age = mean(age, na.rm = TRUE) - qt(0.975, df = n() - 1) * sd(age, na.rm = TRUE) / sqrt(n()),
      ci_upper_age = mean(age, na.rm = TRUE) + qt(0.975, df = n() - 1) * sd(age, na.rm = TRUE) / sqrt(n()),
      n = n(), .groups = "drop") %>%
    group_by(transition) %>%
    mutate(
      proportion = n / sum(n),
      ci_lower_prop = (2 * n + 1.96^2 - 1.96 * sqrt(1.96^2 + 4 * n * (1 - proportion))) / (2 * (sum(n) + 1.96^2)),
      ci_upper_prop = (2 * n + 1.96^2 + 1.96 * sqrt(1.96^2 + 4 * n * (1 - proportion))) / (2 * (sum(n) + 1.96^2))
    ) %>%
    ungroup() %>%
    mutate(group = "Treated")
  
  # Combine, format, and display
  summary_table <- bind_rows(untreated_summary, treated_summary) %>%
    arrange(transition, trans_value, group) %>%
    mutate(
      `Proportion (95% CI)`     = fmt(proportion, ci_lower_prop, ci_upper_prop, digits = 3),
      `Mean Age (95% CI)`       = fmt(mean_age, ci_lower_age, ci_upper_age, digits = 1),
      `Mean Life Exp (95% CI)`  = fmt(mean_life_exp, ci_lower_life_exp, ci_upper_life_exp, digits = 2)
    ) %>%
    select(Group = group, Transition = transition, `Transition Value` = trans_value, n,
           `Proportion (95% CI)`, `Mean Age (95% CI)`, `Mean Life Exp (95% CI)`)
  
  print(summary_table)
  
  write.csv(summary_table, "life_exp_summary_table.csv", row.names = FALSE)
  
  
  # Exit age
  
  library(dplyr)
  library(tidyr)
  
  # Lookup table
  transition_lookup <- c(
    "5"  = "K -> KD",
    "6"  = "K -> KH",
    "7"  = "K -> Death",
    "8"  = "D -> DK",
    "9"  = "D -> DH",
    "10" = "D -> Death",
    "11" = "H -> HK",
    "12" = "H -> HD",
    "13" = "H -> Death",
    "14" = "KD -> KDH",
    "15" = "KD -> Death",
    "16" = "KH -> KHD",
    "17" = "KH -> Death",
    "18" = "DK -> DKH",
    "19" = "DK -> Death",
    "20" = "DH -> DHK",
    "21" = "DH -> Death",
    "22" = "HK -> HKD",
    "23" = "HK -> Death",
    "24" = "HD -> HDK",
    "25" = "HD -> Death",
    "26" = "KDH -> Death",
    "27" = "KHD -> Death",
    "28" = "DKH -> Death",
    "29" = "DHK -> Death",
    "30" = "HKD -> Death",
    "31" = "HDK -> Death"
  )
  
  # Add exit age columns
  pop_b <- pop_b %>%
    mutate(
      exit_age2_soc = age2 + time2,
      exit_age2_trt = age2 + time2_trt,
      exit_age3_soc = age3 + time3,
      exit_age3_trt = age3 + time3_trt,
      exit_age4_soc = age4 + time4,
      exit_age4_trt = age4 + time4_trt
    )
  
  # Pivot to long format
  long_df <- bind_rows(
    # trans2
    pop_b %>%
      filter(!is.na(trans2)) %>%
      select(trans2, exit_age2_soc, exit_age2_trt) %>%
      pivot_longer(
        cols = c(exit_age2_soc, exit_age2_trt),
        names_to = "arm",
        values_to = "exit_age"
      ) %>%
      mutate(
        transition_num = trans2,
        transition = transition_lookup[as.character(trans2)],
        arm = ifelse(arm == "exit_age2_soc", "SOC", "SGLT-2 Inhibitor")
      ) %>%
      select(transition_num, transition, arm, exit_age),
    
    # trans3
    pop_b %>%
      filter(!is.na(trans3)) %>%
      select(trans3, exit_age3_soc, exit_age3_trt) %>%
      pivot_longer(
        cols = c(exit_age3_soc, exit_age3_trt),
        names_to = "arm",
        values_to = "exit_age"
      ) %>%
      mutate(
        transition_num = trans3,
        transition = transition_lookup[as.character(trans3)],
        arm = ifelse(arm == "exit_age3_soc", "SOC", "SGLT-2 Inhibitor")
      ) %>%
      select(transition_num, transition, arm, exit_age),
    
    # trans4
    pop_b %>%
      filter(!is.na(trans4)) %>%
      select(trans4, exit_age4_soc, exit_age4_trt) %>%
      pivot_longer(
        cols = c(exit_age4_soc, exit_age4_trt),
        names_to = "arm",
        values_to = "exit_age"
      ) %>%
      mutate(
        transition_num = trans4,
        transition = transition_lookup[as.character(trans4)],
        arm = ifelse(arm == "exit_age4_soc", "SOC", "SGLT-2 Inhibitor")
      ) %>%
      select(transition_num, transition, arm, exit_age)
  )
  
  # Summarise mean exit age and 95% CI by transition and arm
  summary_df <- long_df %>%
    filter(!is.na(exit_age), !is.na(transition)) %>%
    group_by(transition_num, transition, arm) %>%
    summarise(
      n       = n(),
      mean    = mean(exit_age),
      se      = sd(exit_age) / sqrt(n),
      ci_low  = mean - 1.96 * se,
      ci_high = mean + 1.96 * se,
      .groups = "drop"
    ) %>%
    mutate(
      result = sprintf("%.2f (%.2f, %.2f)", mean, ci_low, ci_high)
    ) %>%
    arrange(transition_num, arm)
  
  print(summary_df)
  
  write.csv(summary_df, "mean_age_exit.csv", row.names = FALSE)
  
  # --------------------------
  # FINAL OUTCOMES
  # --------------------------
  
#  results[[b]] <- data.frame(
#    mean_life_exp = mean(pop_b$life.exp_trt),
#    # Censoring at 100 years old
#    rmst = mean(pmin(pop_b$life.exp_trt, 100)) 
#  )
  
  trans_levels2 <- c(5,6,7,8,9,10,11,12,13)
  tab2 <- table(factor(pop_b$trans2_trt, levels = trans_levels2))
  prop2 <- as.numeric(prop.table(tab2))
  
  trans_levels3 <- c(14,15,16,17,18,19,20,21,22,23,24,25)
  tab3 <- table(factor(pop_b$trans3_trt, levels = trans_levels3))
  prop3 <- as.numeric(prop.table(tab3))
  
  trans_levels4 <- c(26,27,28,29,30,31)
  tab4 <- table(factor(pop_b$trans4_trt, levels = trans_levels4))
  prop4 <- as.numeric(prop.table(tab4))
  
  results_trt[[b]] <- data.frame(
    mean_life_exp = mean(pop_b$life.exp_trt),
    rmst = mean(pmin(pop_b$life.exp_trt, 100)),
    mean_t5 = mean(pop_b$time2_trt[pop_b$trans2_trt==5], na.rm=TRUE),
    mean_t6 = mean(pop_b$time2_trt[pop_b$trans2_trt==6], na.rm=TRUE),
    mean_t7 = mean(pop_b$time2_trt[pop_b$trans2_trt==7], na.rm=TRUE),
    mean_t8 = mean(pop_b$time2_trt[pop_b$trans2_trt==8], na.rm=TRUE),
    mean_t9 = mean(pop_b$time2_trt[pop_b$trans2_trt==9], na.rm=TRUE),
    mean_t10 = mean(pop_b$time2_trt[pop_b$trans2_trt==10], na.rm=TRUE),
    mean_t11 = mean(pop_b$time2_trt[pop_b$trans2_trt==11], na.rm=TRUE),
    mean_t12 = mean(pop_b$time2_trt[pop_b$trans2_trt==12], na.rm=TRUE),
    mean_t13 = mean(pop_b$time2_trt[pop_b$trans2_trt==13], na.rm=TRUE),
    mean_t14 = mean(pop_b$time3_trt[pop_b$trans3_trt==14], na.rm=TRUE),
    mean_t15 = mean(pop_b$time3_trt[pop_b$trans3_trt==15], na.rm=TRUE),
    mean_t16 = mean(pop_b$time3_trt[pop_b$trans3_trt==16], na.rm=TRUE),
    mean_t17 = mean(pop_b$time3_trt[pop_b$trans3_trt==17], na.rm=TRUE),
    mean_t18 = mean(pop_b$time3_trt[pop_b$trans3_trt==18], na.rm=TRUE),
    mean_t19 = mean(pop_b$time3_trt[pop_b$trans3_trt==19], na.rm=TRUE),
    mean_t20 = mean(pop_b$time3_trt[pop_b$trans3_trt==20], na.rm=TRUE),
    mean_t21 = mean(pop_b$time3_trt[pop_b$trans3_trt==21], na.rm=TRUE),
    mean_t22 = mean(pop_b$time3_trt[pop_b$trans3_trt==22], na.rm=TRUE),
    mean_t23 = mean(pop_b$time3_trt[pop_b$trans3_trt==23], na.rm=TRUE),
    mean_t24 = mean(pop_b$time3_trt[pop_b$trans3_trt==24], na.rm=TRUE),
    mean_t25 = mean(pop_b$time3_trt[pop_b$trans3_trt==25], na.rm=TRUE),
    mean_t26 = mean(pop_b$time4_trt[pop_b$trans4_trt==26], na.rm=TRUE),
    mean_t27 = mean(pop_b$time4_trt[pop_b$trans4_trt==27], na.rm=TRUE),
    mean_t28 = mean(pop_b$time4_trt[pop_b$trans4_trt==28], na.rm=TRUE),
    mean_t29 = mean(pop_b$time4_trt[pop_b$trans4_trt==29], na.rm=TRUE),
    mean_t30 = mean(pop_b$time4_trt[pop_b$trans4_trt==30], na.rm=TRUE),
    mean_t31 = mean(pop_b$time4_trt[pop_b$trans4_trt==31], na.rm=TRUE),
    t5 = prop2[1],
    t6 = prop2[2],
    t7 = prop2[3],
    t8 = prop2[4],
    t9 = prop2[5],
    t10 = prop2[6],
    t11 = prop2[7],
    t12 = prop2[8],
    t13 = prop2[9], 
    t14 = prop3[1],
    t15 = prop3[2],
    t16 = prop3[3],
    t17 = prop3[4],
    t18 = prop3[5], 
    t19 = prop3[6],
    t20 = prop3[7],
    t21 = prop3[8], 
    t22 = prop3[9],
    t23 = prop3[10],
    t24 = prop3[11],
    t25 = prop3[12], 
    t26 = prop4[1],
    t27 = prop4[2],
    t28 = prop4[3],
    t29 = prop4[4],
    t30 = prop4[5], 
    t31 = prop4[6]
  )
  
  mc_pop_trt[[b]] <- pop_b[, c("time2_trt", "trans2_trt", "time3_trt", "trans3_trt", "time4_trt", "trans4_trt")]
}

res_trt <- do.call(rbind, results_trt)

summary(res_trt)
apply(res_trt, 2, quantile, c(0.025, 0.5, 0.975))


# Difference in RMST

diff <- res_trt$rmst - res$rmst

diff_df <- data.frame(rmst_diff = diff)

# Scaling this for those with the conditions only, algebraically
# diff_cond = diff_uncond / p1
# As p1 is exactly the same in both arms
# From Table 1, the number of individuals transitioning from ICF to K, D, or H is 179,427 + 74,002 + 80,821 = 334,250 out of 1,000,000, so p1 = 0.33425.
# p1 <- 0.33425
# diff_cond <- diff_life_exp / p1
# quantile(diff_cond, c(0.025, 0.5, 0.975), na.rm = TRUE)
# 0.60 years, 95% CI: −0.09, 1.20

library(ggplot2)

p<- ggplot(diff_df, aes(x = rmst_diff)) +
  geom_histogram(bins = 30) +
#  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    title = "MC distribution of RMST difference",
    x = "RMST (treated - untreated) years",
    y = "Count"
  )

#ggsave(
#  filename = paste0("/Users/rkowen/Library/CloudStorage/OneDrive-SwanseaUniversity/Research_Projects/NIHR_REMIT/HTA_Sim/MonteCarlo/RMST_diff.svg"),
#  plot = p,
#  height = 10,
#  width = 10
#)

p

# Life expectancy
quantile(diff_df$rmst_diff, c(0.025, 0.5, 0.975), na.rm = TRUE)


results_trans <- list()
results_meantrans<- list()

for (i in 5:31) {
  
  var <- paste0("t", i)
  var1 <- paste0("mean_t",i)
  
  diff <- res_trt[[var]] - res[[var]]
  diff_mean <- res_trt[[var1]] - res[[var1]]
  
  results_trans[[var]] <- quantile(diff,
                             probs = c(0.025, 0.5, 0.975),
                             na.rm = TRUE)
  
  results_meantrans[[var1]] <- quantile(diff_mean,
                                   probs = c(0.025, 0.5, 0.975),
                                   na.rm = TRUE)
}

results_trans <- do.call(rbind, results_trans)
results_meantrans <- do.call(rbind, results_meantrans)

# Plots

# Trans2

# Mean KM curve + uncertainty
time_grid <- seq(0, 20, by = 1)

interp_list <- list()

for (i in seq_along(mc_pop)) {
  
  pop_u <- mc_pop[[i]]
  pop_t <- mc_pop_trt[[i]]
  
  for (tr in 5:13) {
    
    dat_u <- subset(pop_u, trans2 == tr)
    dat_t <- subset(pop_t, trans2_trt == tr)
    
    if (nrow(dat_u) < 20 || nrow(dat_t) < 20)
      next
    
    sf_u <- survfit(Surv(time2, rep(1, nrow(dat_u))) ~ 1, data = dat_u)
    sf_t <- survfit(Surv(time2_trt, rep(1, nrow(dat_t))) ~ 1, data = dat_t)
    
    # interpolation
    s_u <- approx(sf_u$time, sf_u$surv, xout = time_grid, rule = 2)$y
    s_t <- approx(sf_t$time, sf_t$surv, xout = time_grid, rule = 2)$y
    
    interp_list[[length(interp_list) + 1]] <- data.frame(
      time = time_grid,
      surv = s_u,
      trans = tr,
      iter = i,
      group = "Untreated"
    )
    
    interp_list[[length(interp_list) + 1]] <- data.frame(
      time = time_grid,
      surv = s_t,
      trans = tr,
      iter = i,
      group = "Treated"
    )
  }
}

km_interp <- do.call(rbind, interp_list)

library(dplyr)

km_summary <- km_interp %>%
  group_by(time, trans, group) %>%
  summarise(
    mean_surv = mean(surv, na.rm = TRUE),
    lcl = quantile(surv, 0.025, na.rm = TRUE),
    ucl = quantile(surv, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

p<- ggplot(km_summary, aes(x = time, y = mean_surv, colour = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = group), alpha = 0.2, colour = NA) +
  facet_wrap(~ trans) +
  labs(
    title = "Kaplan-Meier with Monte Carlo Uncertainty",
    y = "Survival probability"
  )

ggsave(
  filename = paste0("/Users/rkowen/Library/CloudStorage/OneDrive-SwanseaUniversity/Research_Projects/NIHR_REMIT/HTA_Sim/MonteCarlo/KM_trans2.svg"),
  plot = p,
  height = 10,
  width = 20
)

# Trans 3

# Mean KM curve + uncertainty
time_grid <- seq(0, 20, by = 1)

interp_list <- list()

for (i in seq_along(mc_pop)) {
  
  pop_u <- mc_pop[[i]]
  pop_t <- mc_pop_trt[[i]]
  
  for (tr in 14:25) {
    
    dat_u <- subset(pop_u, trans3 == tr)
    dat_t <- subset(pop_t, trans3_trt == tr)
    
    if (nrow(dat_u) < 20 || nrow(dat_t) < 20)
      next
    
    sf_u <- survfit(Surv(time3, rep(1, nrow(dat_u))) ~ 1, data = dat_u)
    sf_t <- survfit(Surv(time3_trt, rep(1, nrow(dat_t))) ~ 1, data = dat_t)
    
    # interpolation
    s_u <- approx(sf_u$time, sf_u$surv, xout = time_grid, rule = 2)$y
    s_t <- approx(sf_t$time, sf_t$surv, xout = time_grid, rule = 2)$y
    
    interp_list[[length(interp_list) + 1]] <- data.frame(
      time = time_grid,
      surv = s_u,
      trans = tr,
      iter = i,
      group = "Untreated"
    )
    
    interp_list[[length(interp_list) + 1]] <- data.frame(
      time = time_grid,
      surv = s_t,
      trans = tr,
      iter = i,
      group = "Treated"
    )
  }
}

km_interp <- do.call(rbind, interp_list)

library(dplyr)

km_summary <- km_interp %>%
  group_by(time, trans, group) %>%
  summarise(
    mean_surv = mean(surv, na.rm = TRUE),
    lcl = quantile(surv, 0.025, na.rm = TRUE),
    ucl = quantile(surv, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

p<- ggplot(km_summary, aes(x = time, y = mean_surv, colour = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = group), alpha = 0.2, colour = NA) +
  facet_wrap(~ trans) +
  labs(
    title = "Kaplan-Meier with Monte Carlo Uncertainty",
    y = "Survival probability"
  )

ggsave(
  filename = paste0("/Users/rkowen/Library/CloudStorage/OneDrive-SwanseaUniversity/Research_Projects/NIHR_REMIT/HTA_Sim/MonteCarlo/KM_trans3.svg"),
  plot = p,
  height = 10,
  width = 20
)

# Trans 4

# Mean KM curve + uncertainty
time_grid <- seq(0, 20, by = 1)

interp_list <- list()

for (i in seq_along(mc_pop)) {
  
  pop_u <- mc_pop[[i]]
  pop_t <- mc_pop_trt[[i]]
  
  for (tr in 26:31) {
    
    dat_u <- subset(pop_u, trans4 == tr)
    dat_t <- subset(pop_t, trans4_trt == tr)
    
    if (nrow(dat_u) < 20 || nrow(dat_t) < 20)
      next
    
    sf_u <- survfit(Surv(time4, rep(1, nrow(dat_u))) ~ 1, data = dat_u)
    sf_t <- survfit(Surv(time4_trt, rep(1, nrow(dat_t))) ~ 1, data = dat_t)
    
    # interpolation
    s_u <- approx(sf_u$time, sf_u$surv, xout = time_grid, rule = 2)$y
    s_t <- approx(sf_t$time, sf_t$surv, xout = time_grid, rule = 2)$y
    
    interp_list[[length(interp_list) + 1]] <- data.frame(
      time = time_grid,
      surv = s_u,
      trans = tr,
      iter = i,
      group = "Untreated"
    )
    
    interp_list[[length(interp_list) + 1]] <- data.frame(
      time = time_grid,
      surv = s_t,
      trans = tr,
      iter = i,
      group = "Treated"
    )
  }
}

km_interp <- do.call(rbind, interp_list)

library(dplyr)

km_summary <- km_interp %>%
  group_by(time, trans, group) %>%
  summarise(
    mean_surv = mean(surv, na.rm = TRUE),
    lcl = quantile(surv, 0.025, na.rm = TRUE),
    ucl = quantile(surv, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

p<- ggplot(km_summary, aes(x = time, y = mean_surv, colour = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = group), alpha = 0.2, colour = NA) +
  facet_wrap(~ trans) +
  labs(
    title = "Kaplan-Meier with Monte Carlo Uncertainty",
    y = "Survival probability"
  )

ggsave(
  filename = paste0("/Users/rkowen/Library/CloudStorage/OneDrive-SwanseaUniversity/Research_Projects/NIHR_REMIT/HTA_Sim/MonteCarlo/KM_trans4.svg"),
  plot = p,
  height = 10,
  width = 20
)
