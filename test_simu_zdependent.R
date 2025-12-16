# ===================================================================
# 0. ENVIRONMENT SETUP AND FILE LOADING
# ===================================================================
# alpha_1, alpha_2,..., alpha_n 
# lambda_ij = O(n^{2-Csparse*2})
rm(list = ls())
library(matrixcalc)
library(Rcpp)
library(RcppEigen)
library(MASS)
library(ggplot2)
library(latex2exp)
library(dplyr)

# source(file = "nonParamteric.R")
# sourceCpp(file = "src/splinecalc.cpp")
# sourceCpp(file = "src/cumcalc.cpp")
source(file = "R/tGenerateC.R") 
source(file = "R/xConstruct.R")
source(file = "codes/analyze_and_plot.R")
# sourceCpp(file = "src/SimSet2.cpp")
sourceCpp(file = "src/SimSet_new.cpp")
# source(file = "sparse_codes2/convert_zij_to_sparse_df.R")
# source(file = "codes/nonParametric_diff_flexible.R")
# sourceCpp("codes/create_tensors_flexible.cpp")
# sourceCpp("codes/compute_flexible.cpp")
source(file = "codes/nonParametric_diff_flexible_2.R")
sourceCpp("codes/create_tensors_flexible_2.cpp")
sourceCpp("codes/compute_flexible_2.cpp")
source(file = "codes/convert_zij_to_sparse_df.R")

set.seed(2025)
n = 100
p = 1
K = 100
rep = 100 # Reduced for quick demonstration
rho = 0.5
sm = T
C_lambda = exp(4); Csparse = 0.5
h_bandwidth = 0.2*(n)^{-0.2} # Kernel smoothing bandwidth
t_grid = seq(1/K, 1, 1/K) %>% round(2)

genzij2<- function(n, rho) {
  # Create the covariance matrix
  Sigma <- matrix(0, n, n)
  for (s in 1:n) {
    for (t in 1:n) {
      if (s == t) {
        Sigma[s, t] <- 3  # Variance
      } else if (abs(s - t) <= 3) {
        Sigma[s, t] <- rho  # Covariance for |s-t| <= 3
      }
    }
  }
  # Generate random vector from multivariate normal distribution
  random_vector <- mvrnorm(1, mu = rep(0, n), Sigma = Sigma)
  return(random_vector)
}
# ===================================================================
# 2. INITIALIZATION AND SIMULATION LOOP
# ===================================================================
# set.seed(123)
# CRITICAL STEP: Generate AND center the covariates
# zij <- array(rbinom(n*n*p, 1, 0.01), c(n, n, p))
# aij <- sample(1:(n*n*p), n*n*p*rho, replace = F)
# zij <- vector(length = n*n*p)
# zij[aij] = runif(n*n*p*rho, -0.5, 0.5)
# zij <- array(zij, c(n, n, p))
zij = array(0, c(n,n,p));
for(i in 1:n){
  # zij[i,,1]=(genzij1(n,rho)>0);
  zij[i,,1]= genzij2(n,rho);
}
ifelse(zij > 0.5, 0.5, ifelse(zij < -0.5, -0.5, zij)) -> zij_
# stopifnot(abs(sum(zij)) < 1e-9)
zij_4d <- array(zij_, c(n, n, p, 1))
z_df <- convert_zij_to_sparse_df(zij)

q_track_alpha <- c(1, 2)
q_track_beta <- c(1, 2)
q_track_gamma <- c(1) # Still works even with a single element vector

# --- MODIFICATION 2: Initialize the results list structure dynamically ---
results <- list()
params <- c("lambda0", "alpha", "beta", "gamma")
types <- c("smooth", "integrated")

# Initialize for lambda0 (which is not indexed by q)
results$lambda0 <- list()
for (type in types) {
  results$lambda0[[type]]$hat <- list()
  results$lambda0[[type]]$var <- list()
}

# Initialize for alpha, beta, gamma, creating sub-lists for each tracked index
initialize_param_results <- function(param_name, q_track_vector) {
  param_list <- list()
  for (type in types) {
    param_list[[type]] <- list(hat = list(), var = list())
    for (q in q_track_vector) {
      q_name <- paste0("q_", q)
      param_list[[type]]$hat[[q_name]] <- list()
      param_list[[type]]$var[[q_name]] <- list()
    }
  }
  return(param_list)
}

results$alpha <- initialize_param_results("alpha", q_track_alpha)
results$beta <- initialize_param_results("beta", q_track_beta)
results$gamma <- initialize_param_results("gamma", q_track_gamma)


cat("Starting simulation...\n")
for (i in 1:rep) {
  t0 <- Sys.time()
  event_data <- tGenerateC(n, C_lambda, Csparse, Zij = zij)
  t1 <- Sys.time()
  np_est <- nonParametric_diff_flexible(event_data$trail,
                                        z_df,
                                        n, p, K,
                                        smooth = sm, h_bandwidth)
  # np_est <- nonParametric_diff_flexible(event_data$trail,
  #                                       zij_4d,
  #                                       n, p, K,
  #                                       smooth = sm, h_bandwidth)
  t3 <- Sys.time()
  
  # Store results for lambda0 (unchanged)
  results$lambda0$integrated$hat[[i]] <- np_est$truncate$integrated$lambda0_hat %>% as.vector()
  results$lambda0$integrated$var[[i]] <- np_est$truncate$integrated$lambda0_var %>% as.vector()
  results$lambda0$smooth$hat[[i]] <- np_est$kernel$pointed$lambda0_hat %>% as.vector()
  results$lambda0$smooth$var[[i]] <- np_est$kernel$pointed$lambda0_var %>% as.vector()
  
  # --- MODIFICATION 3: Loop through q_track vectors to store results ---
  
  # Store results for alpha for each tracked index
  for (q in q_track_alpha) {
    q_name <- paste0("q_", q)
    results$alpha$integrated$hat[[q_name]][[i]] <- np_est$truncate$integrated$alpha_hat[q, ] %>% as.vector()
    results$alpha$integrated$var[[q_name]][[i]] <- np_est$truncate$integrated$alpha_var[q, ] %>% as.vector()
    results$alpha$smooth$hat[[q_name]][[i]] <- np_est$kernel$pointed$alpha_hat[q, ] %>% as.vector()
    results$alpha$smooth$var[[q_name]][[i]] <- np_est$kernel$pointed$alpha_var[q, ] %>% as.vector()
  }
  
  # Store results for beta for each tracked index
  for (q in q_track_beta) {
    q_name <- paste0("q_", q)
    results$beta$integrated$hat[[q_name]][[i]] <- np_est$truncate$integrated$beta_hat[q, ] %>% as.vector()
    results$beta$integrated$var[[q_name]][[i]] <- np_est$truncate$integrated$beta_var[q, ] %>% as.vector()
    results$beta$smooth$hat[[q_name]][[i]] <- np_est$kernel$pointed$beta_hat[q, ] %>% as.vector()
    results$beta$smooth$var[[q_name]][[i]] <- np_est$kernel$pointed$beta_var[q, ] %>% as.vector()
  }
  
  # Store results for gamma for each tracked index
  for (q in q_track_gamma) {
    q_name <- paste0("q_", q)
    results$gamma$integrated$hat[[q_name]][[i]] <- np_est$truncate$integrated$gamma_hat[q, ] %>% as.vector()
    results$gamma$integrated$var[[q_name]][[i]] <- np_est$truncate$integrated$gamma_var[q, q, ] %>% as.vector()
    results$gamma$smooth$hat[[q_name]][[i]] <- np_est$kernel$pointed$gamma_hat[q, ] %>% as.vector()
    results$gamma$smooth$var[[q_name]][[i]] <- np_est$kernel$pointed$gamma_var[q, q, ] %>% as.vector()
  }
  
  cat(".")
  if (i %% 10 == 0) cat("Finished repetition", i, "of", rep, "\n")
}
cat("Simulation complete!\n")
# save(results, file = "result/test_simu_1.rdata")

# ===================================================================
# 3. ANALYSIS AND VISUALIZATION
# ===================================================================
# --- Define the specific time points for coverage calculation ---
coverage_points <- c(0.3, 0.5, 0.7) %>% round(2)


# --- Run analysis for each parameter ---

# 1. Lambda0 / Lambda_0 (Baseline)
analyze_and_plot(results$lambda0$smooth$hat, results$lambda0$smooth$var, unlist(lapply(t_grid, bs_new, n=n, C_lambda = C_lambda, Csparse = Csparse)), t_grid, coverage_points,
                 TeX("Smoothed Baseline Intensity $\\lambda_0(t)$"), "Intensity")

analyze_and_plot(results$lambda0$integrated$hat, results$lambda0$integrated$var, bsI_cpp(t_grid, n, C_lambda, Csparse), t_grid, coverage_points,
                 TeX("Integrated Baseline Intensity $\\Lambda_0(t)$"), "Cumulative Intensity")

# 2. Alpha / A (Sender Effect)
for (q in q_track_alpha) {
  q_name <- paste0("q_", q)
  analyze_and_plot(results$alpha$smooth$hat[[q_name]], results$alpha$smooth$var[[q_name]], unlist(lapply(t_grid, fs_new, i = q, n=n, Csparse = Csparse)), t_grid, coverage_points,
                   TeX(paste0("Smoothed Sender Effect $\\alpha_{", q, "}(t)$")), "Intensity")
  analyze_and_plot(results$alpha$integrated$hat[[q_name]], results$alpha$integrated$var[[q_name]], fsI_cpp(t_grid, i = q, n, Csparse), t_grid, coverage_points,
                   TeX(paste0("Integrated Sender Effect $A_{", q, "}(t)$")), "Cumulative Effect")
}
# 3. Beta / B (Receiver Effect)
for (q in q_track_beta) {
  q_name <- paste0("q_", q)
  analyze_and_plot(results$beta$smooth$hat[[q_name]], results$beta$smooth$var[[q_name]], unlist(lapply(t_grid, fr_new, j = q, n=n, Csparse = Csparse)), t_grid, coverage_points,
                   TeX(paste0("Smoothed Receiver Effect $\\beta_{", q, "}(t)$")), "Intensity")
  analyze_and_plot(results$beta$integrated$hat[[q_name]], results$beta$integrated$var[[q_name]], frI_cpp(t_grid, j = q, n, Csparse), t_grid, coverage_points,
                   TeX(paste0("Integrated Receiver Effect $B_{", q, "}(t)$")), "Cumulative Effect")
}
# 4. Gamma / G (Covariate Effect)
for (q in q_track_gamma) {
  q_name <- paste0("q_", q)
  analyze_and_plot(results$gamma$smooth$hat[[q_name]], results$gamma$smooth$var[[q_name]], unlist(lapply(t_grid, fg_new, n=n, Csparse = Csparse)), t_grid, coverage_points,
                   TeX(paste0("Smoothed Covariate Effect $\\gamma_{", q, "}(t)$")), "Effect Size")
  analyze_and_plot(results$gamma$integrated$hat[[q_name]], results$gamma$integrated$var[[q_name]], fgI_cpp(t_grid, n, Csparse), t_grid, coverage_points,
                   TeX(paste0("Integrated Covariate Effect $\\Gamma_{", q, "}(t)$")), "Cumulative Effect")
}

es <- np_est$kernel$pointed$alpha_hat[1, ]
es_sd <- sqrt(np_est$kernel$pointed$alpha_var[1,])

es <- np_est$truncate$integrated$gamma_hat[1, ]
es_sd <- sqrt(np_est$truncate$integrated$gamma_var[1,1,])

es <- np_est$kernel$pointed$gamma_hat[1, ]
es_sd <- sqrt(np_est$kernel$pointed$gamma_var[1,1,])
# 
# es <- np_est$truncate$integrated$alpha_hat[1, ]
# es_sd <- sqrt(np_est$truncate$integrated$alpha_var[1, ] %>% as.vector())

floordf = data.frame(x = np_est$time_points,
                     y = es,
                     yl = es - 1.96 * es_sd,
                     yu = es + 1.96 * es_sd)
key_dates = c(0.39, 0.68, 0.86)
ggplot(floordf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  geom_vline(
    xintercept = key_dates, 
    linetype = "dotted", 
    color = "black", 
    alpha = 0.7
  )
