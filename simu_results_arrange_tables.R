# source("test_simu_1.R")
# source("test_simu_2.R")
# source("test_simu_4.R")
# source("test_simu_5.R")
i = 2

# n = 200, 500, 1000;
n = 500
q_track_alpha <- c(1, 2)
q_track_beta <- c(1, 2)
q_track_gamma <- c(1) # Sti
sourceCpp(file = paste0("src/SimSet", i, ".cpp"))
load(file = paste0("result/test_simu_", i, "/n=", n, ".rdata"))
# --- Helper function for analysis and plotting (NO CHANGES NEEDED HERE) ---
analyze_and_plot <- function(results_hat, results_var, true_func, t_grid, coverage_points, title, ylab) {
  
  hat_matrix <- do.call(rbind, results_hat)
  var_matrix <- do.call(rbind, results_var)
  
  mean(((hat_matrix - matrix(true_func, 
                             nrow = nrow(hat_matrix),
                             ncol = length(t_grid),
                             byrow = T))^2)) -> MISE
  
  cat(paste("--- MISE for:", title, "is: ", MISE, "---\n"))
  
  
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
  
  plot(p)
}

# ===================================================================
# 3. ANALYSIS AND VISUALIZATION
# ===================================================================
# --- Define the specific time points for coverage calculation ---
coverage_points <- c(0.3, 0.5, 0.7) %>% round(2)
K = 100
t_grid = seq(1/K, 1, 1/K) %>% round(2)
# --- Run analysis for each parameter ---

# # 1. Lambda0 / Lambda_0 (Baseline) - Unchanged
# analyze_and_plot(results$lambda0$smooth$hat, results$lambda0$smooth$var, unlist(lapply(t_grid, bs, n=n, Clambda = Clambda)), t_grid, coverage_points,
#                  TeX("Smoothed Baseline Intensity $\\lambda_0(t)$"), "Intensity")
# analyze_and_plot(results$lambda0$integrated$hat, results$lambda0$integrated$var, bsI_cpp(t_grid, n, Clambda = Clambda), t_grid, coverage_points,
#                  TeX("Integrated Baseline Intensity $\\Lambda_0(t)$"), "Cumulative Intensity")

# --- MODIFICATION 4: Loop through tracked indices for plotting ---

# 2. Alpha / A (Sender Effect)
for (q in q_track_alpha) {
  q_name <- paste0("q_", q)
  pdf(file = paste0("result/test_simu_", i, "/alpha_", q, "_" , n , ".pdf"))
  analyze_and_plot(results$alpha$smooth$hat[[q_name]], results$alpha$smooth$var[[q_name]], unlist(lapply(t_grid, fs, i = q, n=n)), t_grid, coverage_points,
                   TeX(paste0("Smoothed Sender Effect $\\alpha_{", q, "}(t)$")), "Intensity")
  dev.off()
  pdf(file = paste0("result/test_simu_", i, "/Alpha_", q, "_" , n , ".pdf"))
  analyze_and_plot(results$alpha$integrated$hat[[q_name]], results$alpha$integrated$var[[q_name]], fsI_cpp(t_grid, i = q, n), t_grid, coverage_points,
                   TeX(paste0("Integrated Sender Effect $A_{", q, "}(t)$")), "Cumulative Effect")
  dev.off()
}

# 3. Beta / B (Receiver Effect)
for (q in q_track_beta) {
  q_name <- paste0("q_", q)
  pdf(file = paste0("result/test_simu_", i, "/beta_", q, "_" , n , ".pdf"))
  analyze_and_plot(results$beta$smooth$hat[[q_name]], results$beta$smooth$var[[q_name]], unlist(lapply(t_grid, fr, j = q, n=n)), t_grid, coverage_points,
                   TeX(paste0("Smoothed Receiver Effect $\\beta_{", q, "}(t)$")), "Intensity")
  dev.off()
  pdf(file = paste0("result/test_simu_", i, "/Beta_", q, "_" , n , ".pdf"))
  analyze_and_plot(results$beta$integrated$hat[[q_name]], results$beta$integrated$var[[q_name]], frI_cpp(t_grid, j = q, n), t_grid, coverage_points,
                   TeX(paste0("Integrated Receiver Effect $B_{", q, "}(t)$")), "Cumulative Effect")
  dev.off()
}

# 4. Gamma / G (Covariate Effect)
for (q in q_track_gamma) {
  q_name <- paste0("q_", q)
  pdf(file = paste0("result/test_simu_", i, "/gamma_", q, "_" , n , ".pdf"))
  analyze_and_plot(results$gamma$smooth$hat[[q_name]], results$gamma$smooth$var[[q_name]], unlist(lapply(t_grid, fg, n=n)), t_grid, coverage_points,
                   TeX(paste0("Smoothed Covariate Effect $\\gamma_{", q, "}(t)$")), "Effect Size")
  dev.off()
  pdf(file = paste0("result/test_simu_", i, "/Gamma_", q, "_" , n , ".pdf"))
  analyze_and_plot(results$gamma$integrated$hat[[q_name]], results$gamma$integrated$var[[q_name]], fgI_cpp(t_grid, n), t_grid, coverage_points,
                   TeX(paste0("Integrated Covariate Effect $\\Gamma_{", q, "}(t)$")), "Cumulative Effect")
  dev.off()
  }

# The rest of your code for specific manual plotting remains unchanged
# ...
