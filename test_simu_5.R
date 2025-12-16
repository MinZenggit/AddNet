# ===================================================================
# 0. ENVIRONMENT SETUP AND FILE LOADING
# // =============================================================================
# // SIMULATION SETTING: All terms are of order O(log(n)/n^2)

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
sourceCpp(file = "src/SimSet5.cpp")
# source(file = "sparse_codes2/convert_zij_to_sparse_df.R")
# source(file = "codes/nonParametric_diff_flexible.R")
# sourceCpp("codes/create_tensors_flexible.cpp")
# sourceCpp("codes/compute_flexible.cpp")
source(file = "codes/nonParametric_diff_flexible_2.R")
sourceCpp("codes/create_tensors_flexible_2.cpp")
sourceCpp("codes/compute_flexible_2.cpp")
source(file = "codes/convert_zij_to_sparse_df.R")

p = 1
K = 100
rep = 100 # Reduced for quick demonstration
rho = 0.3^2
sm = T
Clambda = exp(5); 
times <- c()
for(n in c(100, 200, 500, 1000, 2000)){
  t_grid = seq(1/K, 1, 1/K) %>% round(2)
  h_bandwidth = 0.15*n^{-0.2}
  # ===================================================================
  # 2. INITIALIZATION AND SIMULATION LOOP
  # ===================================================================
  # CRITICAL STEP: Generate AND center the covariates
  aij <- sample(1:(n*n*p), n*n*p*rho, replace = F)
  zij <- vector(length = n*n*p)
  zij[aij] = runif(n*n*p*rho, -20, 20)
  # zij[aij] = rnorm(n*n*p*rho, 0, 5)
  zij <- array(zij, c(n, n, p))
  zij_4d <- array(zij, c(n, n, p, 1))
  z_df <- convert_zij_to_sparse_df(zij)
  
  # --- MODIFICATION 1: Define q_track variables as vectors ---
  # You can now specify multiple indices to track for each parameter
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
  
  tt <- 0
  cat("Starting simulation...\n")
  for (i in 1:rep) {
    t0 <- Sys.time()
    event_data <- tGenerateC(n, Clambda, 0, Zij = zij)
    t1 <- Sys.time()
    np_est <- nonParametric_diff_flexible(event_data$trail,
                                          z_df,
                                          n, p, K,
                                          smooth = sm, h_bandwidth)
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
    tt <- tt + as.numeric(t3-t1)
    cat(".")
    if (i %% 10 == 0) cat("Finished repetition", i, "of", rep, "\n")
  }
  times = c(as.numeric(tt)/rep, times)
  cat("Simulation of n = ", n ,"complete!. Average time:", as.numeric(tt)/rep, "\n")
  # save(results, file = paste0("result/test_simu_5/n=", n, ".rdata"))
}
cat("\n", times)




# --- Helper function for analysis and plotting (NO CHANGES NEEDED HERE) ---
analyze_and_plot <- function(results_hat, results_var, true_func, t_grid, coverage_points, title, ylab) {
  
  hat_matrix <- do.call(rbind, results_hat)
  var_matrix <- do.call(rbind, results_var)
  
  # 1. Calculate Coverage
  coverage_indices <- which(t_grid %in% coverage_points)
  coverage_at_points <- numeric(length(coverage_indices))
  
  if (length(coverage_indices) > 0) {
    for (j in 1:length(coverage_indices)) {
      idx <- coverage_indices[j]
      coverage_at_points[j] <- mean(abs(hat_matrix[, idx] - true_func[idx]) / sqrt(var_matrix[, idx]) < 1.96, na.rm = TRUE)
    }
    cat(paste("--- Coverage Rate Analysis for:", title, "---\n"))
    coverage_df <- data.frame(Time = t_grid[coverage_indices], Coverage = round(coverage_at_points, 3))
    print(coverage_df)
    cat(paste("Mean coverage over these points:", round(mean(coverage_at_points), 3), "\n\n"))
  } else {
    cat(paste("--- No coverage points found in t_grid for:", title, "---\n\n"))
  }
  
  # 2. Prepare data for ggplot
  plot_df <- data.frame(
    t = t_grid,
    mean_est = colMeans(hat_matrix, na.rm = TRUE),
    lower_ci = apply(hat_matrix, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper_ci = apply(hat_matrix, 2, quantile, probs = 0.975, na.rm = TRUE),
    true_val = true_func
  )
  
  # 3. Create plot
  p <- ggplot(plot_df, aes(x = t)) +
    geom_line(aes(y = mean_est, color = "Estimate")) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "blue") +
    geom_line(aes(y = true_val, color = "True"), linetype = "dashed") +
    scale_color_manual(name = "Legend", values = c("Estimate" = "black", "True" = "red")) +
    labs(title = title, x = "Time", y = ylab) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  print(p)
}

# ===================================================================
# 3. ANALYSIS AND VISUALIZATION
# ===================================================================
# --- Define the specific time points for coverage calculation ---
# coverage_points <- c(0.3, 0.5, 0.7) %>% round(2)
# 
# # --- Run analysis for each parameter ---
# 
# # 1. Lambda0 / Lambda_0 (Baseline) - Unchanged
# analyze_and_plot(results$lambda0$smooth$hat, results$lambda0$smooth$var, unlist(lapply(t_grid, bs, n=n, Clambda = Clambda)), t_grid, coverage_points,
#                  TeX("Smoothed Baseline Intensity $\\lambda_0(t)$"), "Intensity")
# analyze_and_plot(results$lambda0$integrated$hat, results$lambda0$integrated$var, bsI_cpp(t_grid, n, Clambda = Clambda), t_grid, coverage_points,
#                  TeX("Integrated Baseline Intensity $\\Lambda_0(t)$"), "Cumulative Intensity")
# 
# # --- MODIFICATION 4: Loop through tracked indices for plotting ---
# 
# # 2. Alpha / A (Sender Effect)
# for (q in q_track_alpha) {
#   q_name <- paste0("q_", q)
#   analyze_and_plot(results$alpha$smooth$hat[[q_name]], results$alpha$smooth$var[[q_name]], unlist(lapply(t_grid, fs, i = q, n=n)), t_grid, coverage_points,
#                    TeX(paste0("Smoothed Sender Effect $\\alpha_{", q, "}(t)$")), "Intensity")
#   analyze_and_plot(results$alpha$integrated$hat[[q_name]], results$alpha$integrated$var[[q_name]], fsI_cpp(t_grid, i = q, n), t_grid, coverage_points,
#                    TeX(paste0("Integrated Sender Effect $A_{", q, "}(t)$")), "Cumulative Effect")
# }
# 
# # 3. Beta / B (Receiver Effect)
# for (q in q_track_beta) {
#   q_name <- paste0("q_", q)
#   analyze_and_plot(results$beta$smooth$hat[[q_name]], results$beta$smooth$var[[q_name]], unlist(lapply(t_grid, fr, j = q, n=n)), t_grid, coverage_points,
#                    TeX(paste0("Smoothed Receiver Effect $\\beta_{", q, "}(t)$")), "Intensity")
#   analyze_and_plot(results$beta$integrated$hat[[q_name]], results$beta$integrated$var[[q_name]], frI_cpp(t_grid, j = q, n), t_grid, coverage_points,
#                    TeX(paste0("Integrated Receiver Effect $B_{", q, "}(t)$")), "Cumulative Effect")
# }
# 
# # 4. Gamma / G (Covariate Effect)
# for (q in q_track_gamma) {
#   q_name <- paste0("q_", q)
#   analyze_and_plot(results$gamma$smooth$hat[[q_name]], results$gamma$smooth$var[[q_name]], unlist(lapply(t_grid, fg, n=n)), t_grid, coverage_points,
#                    TeX(paste0("Smoothed Covariate Effect $\\gamma_{", q, "}(t)$")), "Effect Size")
#   analyze_and_plot(results$gamma$integrated$hat[[q_name]], results$gamma$integrated$var[[q_name]], fgI_cpp(t_grid, n), t_grid, coverage_points,
#                    TeX(paste0("Integrated Covariate Effect $\\Gamma_{", q, "}(t)$")), "Cumulative Effect")
# }
# 
# # The rest of your code for specific manual plotting remains unchanged
# # ...
