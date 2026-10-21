## ============================================================================
## PSA: each of the B iterations draws a fresh SOC
## population (stages 1-4) and a fresh HR_b,
## so SOC and treated are paired within each iteration,
## and pooling across iterations captures both first-order (patient
## sampling) and second-order (HR parameter) uncertainty.
## ============================================================================

rm(list=ls())
library(dplyr)
library(flexsurv)
library(ggplot2)
library(tidyr)
library(survival)
library(scales)  # for percent_format

setwd("wd")

######################
##Simulate a cohort ##
######################
set.seed(1234)

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

probs <- c(0.20, 0.20, 0.20, 0.20, 0.20)

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


# Monte Carlo wrapper #
B <- 100   # number of simulations
n <- nrow(pop)

results <- vector("list", B)
mc_pop      <- vector("list", B)

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

## ---- Per-iteration stat-extraction function ----
compute_iteration_stats <- function(pop_b) {
  rows <- list()
  add_row <- function(stage, label, arm, n_init, sel, age_v, life_v) {
    n_chg <- sum(sel, na.rm = TRUE)
    rows[[length(rows) + 1]] <<- data.frame(
      stage = stage, label = label, arm = arm,
      n_init = n_init, n_chg = n_chg,
      mean_age = if (n_chg > 0) mean(age_v[sel], na.rm = TRUE) else NA_real_,
      mean_life = if (n_chg > 0) mean(life_v[sel], na.rm = TRUE) else NA_real_
    )
  }
  n_cohort <- nrow(pop_b)
  
  stage1_map <- list(K = 1, D = 2, H = 3, Death = 4)
  for (dest in names(stage1_map)) {
    sel <- pop_b$trans1 == stage1_map[[dest]]
    add_row(1, paste("ICF ->", dest), "SOC", n_cohort, sel,
            pop_b$age + pop_b$time1, pop_b$life.exp)
  }
  
  stage2 <- list(
    K = list(parent = 1, children = list(KD = 5, KH = 6, Death = 7)),
    D = list(parent = 2, children = list(DK = 8, DH = 9, Death = 10)),
    H = list(parent = 3, children = list(HK = 11, HD = 12, Death = 13))
  )
  for (node in names(stage2)) {
    n_init <- sum(pop_b$trans1 == stage2[[node]]$parent, na.rm = TRUE)
    for (dest in names(stage2[[node]]$children)) {
      code <- stage2[[node]]$children[[dest]]
      label <- paste(node, "->", dest)
      for (arm in c("SOC", "Treated")) {
        if (arm == "SOC") { tvec <- pop_b$trans2; exitvec <- pop_b$exit_age2_soc; lifevec <- pop_b$life.exp }
        else               { tvec <- pop_b$trans2_trt; exitvec <- pop_b$exit_age2_trt; lifevec <- pop_b$life.exp_trt }
        sel <- tvec == code; sel[is.na(sel)] <- FALSE
        add_row(2, label, arm, n_init, sel, exitvec, lifevec)
      }
    }
  }
  
  stage3 <- list(
    KD = list(parent = 5,  children = list(KDH = 14, Death = 15)),
    KH = list(parent = 6,  children = list(KHD = 16, Death = 17)),
    DK = list(parent = 8,  children = list(DKH = 18, Death = 19)),
    DH = list(parent = 9,  children = list(DHK = 20, Death = 21)),
    HK = list(parent = 11, children = list(HKD = 22, Death = 23)),
    HD = list(parent = 12, children = list(HDK = 24, Death = 25))
  )
  for (node in names(stage3)) {
    for (dest in names(stage3[[node]]$children)) {
      code <- stage3[[node]]$children[[dest]]
      label <- paste(node, "->", dest)
      for (arm in c("SOC", "Treated")) {
        if (arm == "SOC") { pvec <- pop_b$trans2; tvec <- pop_b$trans3; exitvec <- pop_b$exit_age3_soc; lifevec <- pop_b$life.exp }
        else               { pvec <- pop_b$trans2_trt; tvec <- pop_b$trans3_trt; exitvec <- pop_b$exit_age3_trt; lifevec <- pop_b$life.exp_trt }
        n_init <- sum(pvec == stage3[[node]]$parent, na.rm = TRUE)
        sel <- tvec == code; sel[is.na(sel)] <- FALSE
        add_row(3, label, arm, n_init, sel, exitvec, lifevec)
      }
    }
  }
  
  stage4 <- list(KDH = list(parent = 14, code = 26), KHD = list(parent = 16, code = 27),
                 DKH = list(parent = 18, code = 28), DHK = list(parent = 20, code = 29),
                 HKD = list(parent = 22, code = 30), HDK = list(parent = 24, code = 31))
  for (node in names(stage4)) {
    label <- paste(node, "-> Death")
    code <- stage4[[node]]$code
    for (arm in c("SOC", "Treated")) {
      if (arm == "SOC") { pvec <- pop_b$trans3; tvec <- pop_b$trans4; exitvec <- pop_b$exit_age4_soc; lifevec <- pop_b$life.exp }
      else               { pvec <- pop_b$trans3_trt; tvec <- pop_b$trans4_trt; exitvec <- pop_b$exit_age4_trt; lifevec <- pop_b$life.exp_trt }
      n_init <- sum(pvec == stage4[[node]]$parent, na.rm = TRUE)
      sel <- tvec == code; sel[is.na(sel)] <- FALSE
      add_row(4, label, arm, n_init, sel, exitvec, lifevec)
    }
  }
  do.call(rbind, rows)
}

iter_stats <- vector("list", B)
iter_overall <- vector("list", B)

surv_grid  <- 0:20       # years since entry to the state
surv_trans <- 5:31       # transitions with a treatment effect
iter_surv  <- vector("list", B)

# Transition labels (CKD = K, T2DM = D, HF = H), as in Table 2
trans_labels <- c(
  `5`  = "CKD -> CKD-T2DM",            `6`  = "CKD -> CKD-HF",            `7`  = "CKD -> Death",
  `8`  = "T2DM -> T2DM-CKD",           `9`  = "T2DM -> T2DM-HF",          `10` = "T2DM -> Death",
  `11` = "HF -> HF-CKD",               `12` = "HF -> HF-T2DM",            `13` = "HF -> Death",
  `14` = "CKD-T2DM -> CKD-T2DM-HF",    `15` = "CKD-T2DM -> Death",
  `16` = "CKD-HF -> CKD-HF-T2DM",      `17` = "CKD-HF -> Death",
  `18` = "T2DM-CKD -> T2DM-CKD-HF",    `19` = "T2DM-CKD -> Death",
  `20` = "T2DM-HF -> T2DM-HF-CKD",     `21` = "T2DM-HF -> Death",
  `22` = "HF-CKD -> HF-CKD-T2DM",      `23` = "HF-CKD -> Death",
  `24` = "HF-T2DM -> HF-T2DM-CKD",     `25` = "HF-T2DM -> Death",
  `26` = "CKD-T2DM-HF -> Death",       `27` = "CKD-HF-T2DM -> Death",
  `28` = "T2DM-CKD-HF -> Death",       `29` = "T2DM-HF-CKD -> Death",
  `30` = "HF-CKD-T2DM -> Death",       `31` = "HF-T2DM-CKD -> Death"
)


for (b in 1:B) {
  
  cat("Iteration:", b, "\n")
  
  pop_b <- pop
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
  
  
  pop_b$lp9 <- mu0 + beta_age*pop_b$agec2 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
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
  
  pop_b$lp11 <- mu0 + beta_age*pop_b$agec2 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
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
  
  pop_b$lp12 <- mu0 + beta_age*pop_b$agec2 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
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
  
  pop_b$lp14<- log(rate) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
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
  
  pop_b$lp16<- log(rate) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
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
  
  pop_b$lp22<- log(rate) + beta_age*pop_b$agec3 + beta_sex*pop_b$sex + beta_wimd1*pop_b$wimd1 + beta_wimd2*pop_b$wimd2 + beta_wimd4*pop_b$wimd4 + beta_wimd5*pop_b$wimd5
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
  beta_sex<- 0.054890
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
  
  HR_b <- sample_HR(hr_table)   # fresh draw every iteration - proper PSA
  
  # --------------------------
  # STAGE 1–4
  # --------------------------
  
  HR <- HR_b[1]
  sf <- survfit(Surv(t5, rep(1,n)) ~ 1, data = pop_b)
  u5 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t5)  # patient's own SOC percentile (rank-preserving)
  pop_b$t5_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u5)
  
  HR <- HR_b[2]
  sf <- survfit(Surv(t6, rep(1,n)) ~ 1, data = pop_b)
  u6 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t6)  # patient's own SOC percentile (rank-preserving)
  pop_b$t6_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u6)
  
  HR <- HR_b[3]
  sf <- survfit(Surv(t7, rep(1,n)) ~ 1, data = pop_b)
  u7 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t7)  # patient's own SOC percentile (rank-preserving)
  pop_b$t7_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u7)
  
  HR <- HR_b[4]
  sf <- survfit(Surv(t8, rep(1,n)) ~ 1, data = pop_b)
  u8 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t8)  # patient's own SOC percentile (rank-preserving)
  pop_b$t8_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u8)
  
  HR <- HR_b[5]
  sf <- survfit(Surv(t9, rep(1,n)) ~ 1, data = pop_b)
  u9 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t9)  # patient's own SOC percentile (rank-preserving)
  pop_b$t9_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u9)
  
  HR <- HR_b[6]
  sf <- survfit(Surv(t10, rep(1,n)) ~ 1, data = pop_b)
  u10 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t10)  # patient's own SOC percentile (rank-preserving)
  pop_b$t10_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u10)
  
  HR <- HR_b[7]
  sf <- survfit(Surv(t11, rep(1,n)) ~ 1, data = pop_b)
  u11 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t11)  # patient's own SOC percentile (rank-preserving)
  pop_b$t11_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u11)
  
  HR <- HR_b[8]
  sf <- survfit(Surv(t12, rep(1,n)) ~ 1, data = pop_b)
  u12 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t12)  # patient's own SOC percentile (rank-preserving)
  pop_b$t12_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u12)
  
  HR <- HR_b[9]
  sf <- survfit(Surv(t13, rep(1,n)) ~ 1, data = pop_b)
  u13 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t13)  # patient's own SOC percentile (rank-preserving)
  pop_b$t13_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u13)
  
  HR <- HR_b[10]
  sf <- survfit(Surv(t14, rep(1,n)) ~ 1, data = pop_b)
  u14 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t14)  # patient's own SOC percentile (rank-preserving)
  pop_b$t14_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u14)
  
  HR <- HR_b[11]
  sf <- survfit(Surv(t15, rep(1,n)) ~ 1, data = pop_b)
  u15 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t15)  # patient's own SOC percentile (rank-preserving)
  pop_b$t15_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u15)
  
  HR <- HR_b[12]
  sf <- survfit(Surv(t16, rep(1,n)) ~ 1, data = pop_b)
  u16 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t16)  # patient's own SOC percentile (rank-preserving)
  pop_b$t16_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u16)
  
  HR <- HR_b[13]
  sf <- survfit(Surv(t17, rep(1,n)) ~ 1, data = pop_b)
  u17 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t17)  # patient's own SOC percentile (rank-preserving)
  pop_b$t17_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u17)
  
  HR <- HR_b[14]
  sf <- survfit(Surv(t18, rep(1,n)) ~ 1, data = pop_b)
  u18 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t18)  # patient's own SOC percentile (rank-preserving)
  pop_b$t18_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u18)
  
  HR <- HR_b[15]
  sf <- survfit(Surv(t19, rep(1,n)) ~ 1, data = pop_b)
  u19 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t19)  # patient's own SOC percentile (rank-preserving)
  pop_b$t19_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u19)
  
  HR <- HR_b[16]
  sf <- survfit(Surv(t20, rep(1,n)) ~ 1, data = pop_b)
  u20 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t20)  # patient's own SOC percentile (rank-preserving)
  pop_b$t20_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u20)
  
  HR <- HR_b[17]
  sf <- survfit(Surv(t21, rep(1,n)) ~ 1, data = pop_b)
  u21 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t21)  # patient's own SOC percentile (rank-preserving)
  pop_b$t21_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u21)
  
  HR <- HR_b[18]
  sf <- survfit(Surv(t22, rep(1,n)) ~ 1, data = pop_b)
  u22 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t22)  # patient's own SOC percentile (rank-preserving)
  pop_b$t22_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u22)
  
  HR <- HR_b[19]
  sf <- survfit(Surv(t23, rep(1,n)) ~ 1, data = pop_b)
  u23 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t23)  # patient's own SOC percentile (rank-preserving)
  pop_b$t23_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u23)
  
  HR <- HR_b[20]
  sf <- survfit(Surv(t24, rep(1,n)) ~ 1, data = pop_b)
  u24 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t24)  # patient's own SOC percentile (rank-preserving)
  pop_b$t24_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u24)
  
  HR <- HR_b[21]
  sf <- survfit(Surv(t25, rep(1,n)) ~ 1, data = pop_b)
  u25 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t25)  # patient's own SOC percentile (rank-preserving)
  pop_b$t25_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u25)
  
  HR <- HR_b[22]
  sf <- survfit(Surv(t26, rep(1,n)) ~ 1, data = pop_b)
  u26 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t26)  # patient's own SOC percentile (rank-preserving)
  pop_b$t26_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u26)
  
  HR <- HR_b[23]
  sf <- survfit(Surv(t27, rep(1,n)) ~ 1, data = pop_b)
  u27 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t27)  # patient's own SOC percentile (rank-preserving)
  pop_b$t27_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u27)
  
  HR <- HR_b[24]
  sf <- survfit(Surv(t28, rep(1,n)) ~ 1, data = pop_b)
  u28 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t28)  # patient's own SOC percentile (rank-preserving)
  pop_b$t28_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u28)
  
  HR <- HR_b[25]
  sf <- survfit(Surv(t29, rep(1,n)) ~ 1, data = pop_b)
  u29 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t29)  # patient's own SOC percentile (rank-preserving)
  pop_b$t29_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u29)
  
  HR <- HR_b[26]
  sf <- survfit(Surv(t30, rep(1,n)) ~ 1, data = pop_b)
  u30 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t30)  # patient's own SOC percentile (rank-preserving)
  pop_b$t30_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u30)
  
  HR <- HR_b[27]
  sf <- survfit(Surv(t31, rep(1,n)) ~ 1, data = pop_b)
  u31 <- approxfun(sf$time, sf$surv, rule=2)(pop_b$t31)  # patient's own SOC percentile (rank-preserving)
  pop_b$t31_trt <- approxfun(sf$surv^HR, sf$time, rule=2)(u31)

  iter_surv[[b]] <- do.call(rbind, lapply(surv_trans, function(k) {
    t_soc <- pop_b[[paste0("t", k)]]
    t_trt <- pop_b[[paste0("t", k, "_trt")]]
    data.frame(
      iter  = b,
      trans = k,
      time  = rep(surv_grid, 2),
      arm   = rep(c("SOC", "SGLT-2"), each = length(surv_grid)),
      surv  = c(1 - ecdf(t_soc)(surv_grid),   # S(t) = P(T > t)
                1 - ecdf(t_trt)(surv_grid))
    )
  }))

  
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
  

  pop_b <- pop_b %>%
    mutate(
      time2_trt_age = ifelse(is.na(time2_trt), 0, time2_trt),
      age3_trt      = age2 + time2_trt_age,
      time3_trt_age = ifelse(is.na(time3_trt), 0, time3_trt),
      age4_trt      = age3_trt + time3_trt_age,
      exit_age2_soc = age2 + time2,
      exit_age2_trt = age2 + time2_trt,
      exit_age3_soc = age3 + time3,
      exit_age3_trt = age3_trt + time3_trt,
      exit_age4_soc = age4 + time4,
      exit_age4_trt = age4_trt + time4_trt
    )

    iter_state[[b]] <- do.call(rbind, lapply(names(state_defs), function(s) {
    d <- state_defs[[s]]
    do.call(rbind, lapply(c("SOC", "SGLT-2"), function(arm) {
      sfx <- if (arm == "SOC") "" else "_trt"
      # who entered the state in this arm (stage-1 entry is the same in both arms)
      entry_col <- if (d$stage == 1) "trans1" else paste0("trans", d$stage, sfx)
      in_state  <- which(pop_b[[entry_col]] == d$code)
      if (length(in_state) == 0) return(NULL)
      # competing exit times (time since entry) and which exit happens first
      tm   <- sapply(d$exits, function(k) pop_b[[paste0("t", k, sfx)]][in_state])
      tm   <- matrix(tm, ncol = length(d$exits))
      t_ex <- apply(tm, 1, min)
      first <- d$exits[max.col(-tm, ties.method = "first")]
      rbind(
        data.frame(iter = b, state = s, arm = arm, curve = "Remaining in state",
                   time = surv_grid, value = 1 - ecdf(t_ex)(surv_grid)),
        do.call(rbind, lapply(d$exits, function(k)
          data.frame(iter = b, state = s, arm = arm, curve = trans_labels[as.character(k)],
                     time = surv_grid,
                     value = sapply(surv_grid, function(g) mean(t_ex <= g & first == k)))))
      )
    }))
  }))
        
  iter_stats[[b]] <- compute_iteration_stats(pop_b)
  iter_stats[[b]]$iter <- b
  
  # eligible = developed at least one of T2DM/CKD/HF, i.e. ever exposed to
  # the modelled treatment effect (SGLT-2 inhibitors are only indicated
  # once one of these conditions is present)
  eligible <- pop_b$trans1 %in% c(1, 2, 3)
  
  iter_overall[[b]] <- data.frame(
    iter = b,
    mean_life_soc = mean(pop_b$life.exp),
    mean_life_trt = mean(pop_b$life.exp_trt),
    diff_life     = mean(pop_b$life.exp_trt) - mean(pop_b$life.exp),
    n_eligible             = sum(eligible),
    mean_life_soc_eligible = mean(pop_b$life.exp[eligible]),
    mean_life_trt_eligible = mean(pop_b$life.exp_trt[eligible]),
    diff_life_eligible     = mean(pop_b$life.exp_trt[eligible]) - mean(pop_b$life.exp[eligible])
  )
}

## ---- Pool across all B iterations ----
pooled <- do.call(rbind, iter_stats)
pooled$prob <- pooled$n_chg / pooled$n_init
pooled_overall <- do.call(rbind, iter_overall)

key <- interaction(pooled$stage, pooled$label, pooled$arm, drop = TRUE)
pooled_summary <- do.call(rbind, lapply(split(pooled, key), function(g) {
  data.frame(
    stage = g$stage[1], label = g$label[1], arm = g$arm[1],
    n_iter = nrow(g),
    mean_n_init = mean(g$n_init), mean_n_chg = mean(g$n_chg),
    prob_est = mean(g$prob, na.rm = TRUE),
    prob_lo  = quantile(g$prob, 0.025, na.rm = TRUE, names = FALSE),
    prob_hi  = quantile(g$prob, 0.975, na.rm = TRUE, names = FALSE),
    age_est  = mean(g$mean_age, na.rm = TRUE),
    age_lo   = quantile(g$mean_age, 0.025, na.rm = TRUE, names = FALSE),
    age_hi   = quantile(g$mean_age, 0.975, na.rm = TRUE, names = FALSE),
    life_est = mean(g$mean_life, na.rm = TRUE),
    life_lo  = quantile(g$mean_life, 0.025, na.rm = TRUE, names = FALSE),
    life_hi  = quantile(g$mean_life, 0.975, na.rm = TRUE, names = FALSE)
  )
}))
rownames(pooled_summary) <- NULL

## ---- Pool the OVERALL population-level effect 
overall_est <- mean(pooled_overall$diff_life)
overall_ci  <- quantile(pooled_overall$diff_life, c(0.025, 0.975))
cat(sprintf("Overall pooled life expectancy difference: %.2f years (95%% CI: %.2f, %.2f)\n",
            overall_est, overall_ci[1], overall_ci[2]))

## ---- Pool the treatment-ELIGIBLE population's life expectancy effect
## (restricted to those who develop at least one of T2DM/CKD/HF) 
eligible_est <- mean(pooled_overall$diff_life_eligible)
eligible_ci  <- quantile(pooled_overall$diff_life_eligible, c(0.025, 0.975))
cat(sprintf("Treatment-eligible population life expectancy difference: %.2f years (95%% CI: %.2f, %.2f)\n",
            eligible_est, eligible_ci[1], eligible_ci[2]))
write.csv(pooled_overall, "pooled_overall_incl_eligible.csv", row.names = FALSE)

## ---- Pool the per-transition DIFFERENCES 
soc_wide <- pooled[pooled$arm == "SOC", c("stage","label","iter","prob","mean_age","mean_life")]
trt_wide <- pooled[pooled$arm == "Treated", c("stage","label","iter","prob","mean_age","mean_life")]
diffs <- merge(soc_wide, trt_wide, by = c("stage","label","iter"), suffixes = c("_soc","_trt"))
diffs$prob_diff <- diffs$prob_trt - diffs$prob_soc
diffs$age_diff  <- diffs$mean_age_trt - diffs$mean_age_soc
diffs$life_diff <- diffs$mean_life_trt - diffs$mean_life_soc
write.csv(diffs, "per_iteration_differences.csv", row.names = FALSE)  # per-iteration paired differences, used by MC_stability.R

key2 <- interaction(diffs$stage, diffs$label, drop = TRUE)
diff_summary <- do.call(rbind, lapply(split(diffs, key2), function(g) {
  data.frame(
    stage = g$stage[1], label = g$label[1], n_iter = nrow(g),
    prob_diff_est = mean(g$prob_diff, na.rm=TRUE), prob_diff_lo = quantile(g$prob_diff, 0.025, names = FALSE), prob_diff_hi = quantile(g$prob_diff, 0.975, names = FALSE),
    age_diff_est  = mean(g$age_diff, na.rm=TRUE),  age_diff_lo  = quantile(g$age_diff, 0.025, na.rm=TRUE, names = FALSE),  age_diff_hi  = quantile(g$age_diff, 0.975, na.rm=TRUE, names = FALSE),
    life_diff_est = mean(g$life_diff, na.rm=TRUE), life_diff_lo = quantile(g$life_diff, 0.025, na.rm=TRUE, names = FALSE), life_diff_hi = quantile(g$life_diff, 0.975, na.rm=TRUE, names = FALSE)
  )
}))
rownames(diff_summary) <- NULL

print(pooled_summary, digits = 4)
print(diff_summary, digits = 4)
write.csv(pooled_summary, "pooled_transition_table.csv", row.names = FALSE)
write.csv(diff_summary, "pooled_difference_table.csv", row.names = FALSE)

surv_all <- do.call(rbind, iter_surv)
write.csv(surv_all, "per_iteration_survival_curves.csv", row.names = FALSE)

# Mean curve and 2.5th / 97.5th percentiles across the B iterations
surv_summary <- surv_all %>%
  group_by(trans, arm, time) %>%
  summarise(
    est   = mean(surv),
    lower = quantile(surv, 0.025),
    upper = quantile(surv, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    label = factor(trans_labels[as.character(trans)], levels = trans_labels),
    arm   = factor(arm, levels = c("SGLT-2", "SOC"))
  )

plot_surv <- function(trans_set) {
  ggplot(filter(surv_summary, trans %in% trans_set),
         aes(x = time, y = est, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ label, ncol = 3) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Time (years)", y = "Survival probability",
         colour = NULL, fill = NULL) +
    theme_classic() +
    theme(strip.background = element_rect(colour = "black"),
          legend.position = "right")
}

# fig3A <- plot_surv(5:13)    # A: transitions from single conditions
# fig3B <- plot_surv(14:25)   # B: transitions from two conditions
# fig3C <- plot_surv(26:31)   # C: transitions from three conditions

# ggsave("Figure3A_single_conditions.png", fig3A, width = 10, height = 7.5, dpi = 300)
# ggsave("Figure3B_two_conditions.png",    fig3B, width = 10, height = 10,  dpi = 300)
# ggsave("Figure3C_three_conditions.png",  fig3C, width = 10, height = 5,   dpi = 300)

library(dplyr)
library(ggplot2)

surv_all <- read.csv("per_iteration_survival_curves.csv")   # iter, trans, time, arm, surv

trans_labels <- c(
  `5`  = "CKD -> CKD-T2DM",            `6`  = "CKD -> CKD-HF",            `7`  = "CKD -> Death",
  `8`  = "T2DM -> T2DM-CKD",           `9`  = "T2DM -> T2DM-HF",          `10` = "T2DM -> Death",
  `11` = "HF -> HF-CKD",               `12` = "HF -> HF-T2DM",            `13` = "HF -> Death",
  `14` = "CKD-T2DM -> CKD-T2DM-HF",    `15` = "CKD-T2DM -> Death",
  `16` = "CKD-HF -> CKD-HF-T2DM",      `17` = "CKD-HF -> Death",
  `18` = "T2DM-CKD -> T2DM-CKD-HF",    `19` = "T2DM-CKD -> Death",
  `20` = "T2DM-HF -> T2DM-HF-CKD",     `21` = "T2DM-HF -> Death",
  `22` = "HF-CKD -> HF-CKD-T2DM",      `23` = "HF-CKD -> Death",
  `24` = "HF-T2DM -> HF-T2DM-CKD",     `25` = "HF-T2DM -> Death",
  `26` = "CKD-T2DM-HF -> Death",       `27` = "CKD-HF-T2DM -> Death",
  `28` = "T2DM-CKD-HF -> Death",       `29` = "T2DM-HF-CKD -> Death",
  `30` = "HF-CKD-T2DM -> Death",       `31` = "HF-T2DM-CKD -> Death"
)

# 1 - S(t) in each iteration, then mean and 2.5th / 97.5th percentiles across iterations
# (taking percentiles of 1 - S is the same as 1 - the percentiles of S, with lower/upper swapped)
cum_summary <- surv_all %>%
  mutate(cum = 1 - surv) %>%
  group_by(trans, arm, time) %>%
  summarise(
    est   = mean(cum),
    lower = quantile(cum, 0.025, names = FALSE),
    upper = quantile(cum, 0.975, names = FALSE),
    .groups = "drop"
  ) %>%
  mutate(
    label = factor(unname(trans_labels[as.character(trans)]), levels = unname(trans_labels)),
    arm   = factor(arm, levels = c("SGLT-2", "SOC"))
  )

write.csv(cum_summary, "Figure_one_minus_S_summary.csv", row.names = FALSE)

plot_cum <- function(trans_set) {
  ggplot(filter(cum_summary, trans %in% trans_set),
         aes(x = time, y = est, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ label, ncol = 3) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Time (years)", y = "Cumulative probability, 1 - S(t)",
         colour = NULL, fill = NULL) +
    theme_classic() +
    theme(strip.background = element_rect(colour = "black"),
          legend.position = "right")
}

# fig3A_cum <- plot_cum(5:13)    # A: transitions from single conditions
# fig3B_cum <- plot_cum(14:25)   # B: transitions from two conditions
# fig3C_cum <- plot_cum(26:31)   # C: transitions from three conditions

# ggsave("Figure3A_single_conditions_1minusS.png", fig3A_cum, width = 10, height = 7.5, dpi = 300)
# ggsave("Figure3B_two_conditions_1minusS.png",    fig3B_cum, width = 10, height = 10,  dpi = 300)
# ggsave("Figure3C_three_conditions_1minusS.png",  fig3C_cum, width = 10, height = 5,   dpi = 300)

state_all <- do.call(rbind, iter_state)
write.csv(state_all, "per_iteration_state_curves.csv", row.names = FALSE)

state_summary <- state_all %>%
  group_by(state, arm, curve, time) %>%
  summarise(est = mean(value), lower = quantile(value, 0.025),
            upper = quantile(value, 0.975), .groups = "drop") %>%
  mutate(arm = factor(arm, levels = c("SGLT-2", "SOC")))

# (i) Probability of remaining in each state - one panel per state
plot_remain <- function(states) {
  ggplot(filter(state_summary, curve == "Remaining in state", state %in% states) %>%
           mutate(state = factor(state, levels = states)),
         aes(time, est, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ state, ncol = 3) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Time since entering state (years)",
         y = "Probability of remaining in state", colour = NULL, fill = NULL) +
    theme_classic() + theme(strip.background = element_rect(colour = "black"))
}

# (ii) Cumulative incidence of each exit - one panel per transition
plot_cif <- function(states) {
  ggplot(filter(state_summary, curve != "Remaining in state", state %in% states) %>%
           mutate(curve = factor(curve, levels = trans_labels)),
         aes(time, est, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ curve, ncol = 3) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Time since entering state (years)",
         y = "Cumulative incidence of transition", colour = NULL, fill = NULL) +
    theme_classic() + theme(strip.background = element_rect(colour = "black"))
}

single <- c("CKD", "T2DM", "HF")
double <- c("CKD-T2DM", "CKD-HF", "T2DM-CKD", "T2DM-HF", "HF-CKD", "HF-T2DM")
triple <- c("CKD-T2DM-HF", "CKD-HF-T2DM", "T2DM-CKD-HF", "T2DM-HF-CKD", "HF-CKD-T2DM", "HF-T2DM-CKD")

ggsave("Figure3A_remaining_single.png", plot_remain(single), width = 10, height = 3.5, dpi = 300)
ggsave("Figure3B_remaining_double.png", plot_remain(double), width = 10, height = 6,   dpi = 300)
ggsave("Figure3C_remaining_triple.png", plot_remain(triple), width = 10, height = 6,   dpi = 300)

ggsave("Figure3A_cif_single.png", plot_cif(single), width = 10, height = 7.5, dpi = 300)
ggsave("Figure3B_cif_double.png", plot_cif(double), width = 10, height = 10,  dpi = 300)
ggsave("Figure3C_cif_triple.png", plot_cif(triple), width = 10, height = 5, dpi = 300)
