# --- Helper function for analysis and plotting ---
# *** MODIFIED to use your robust, point-by-point coverage calculation method ***
analyze_and_plot <- function(results_hat, results_var, true_func, t_grid, coverage_points, title, ylab) {
  
  hat_matrix <- do.call(rbind, results_hat)
  var_matrix <- do.call(rbind, results_var)
  
  apply(hat_matrix, 2, mean)
  apply(hat_matrix, 2, var)
  apply(var_matrix, 2, mean)
  # 1. Calculate Coverage using a loop for robustness (as you suggested)
  coverage_indices <- which(t_grid %in% coverage_points)
  coverage_at_points <- numeric(length(coverage_indices))
  
  if (length(coverage_indices) > 0) {
    for (j in 1:length(coverage_indices)) {
      idx <- coverage_indices[j]
      # This is the robust calculation for a single point `t`
      coverage_at_points[j] <- mean(abs(hat_matrix[, idx] - true_func[idx]) / sqrt(var_matrix[, idx]) < 1.96, na.rm = TRUE)
    }
    
    # Report coverage
    cat(paste("--- Coverage Rate Analysis for:", title, "---\n"))
    coverage_df <- data.frame(Time = t_grid[coverage_indices], Coverage = round(coverage_at_points, 3))
    print(coverage_df)
    cat(paste("Mean coverage over these points:", round(mean(coverage_at_points), 3), "\n\n"))
  } else {
    cat(paste("--- No coverage points found in t_grid for:", title, "---\n\n"))
  }
  
  # 2. Prepare data for ggplot (using the full grid)
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
