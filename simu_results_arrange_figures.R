# ===================================================================
# 0. SETUP - 加载必要的库和定义全局变量
# ===================================================================
# install.packages(c("dplyr", "tidyr", "knitr", "ggplot2", "Rcpp"))
library(dplyr)
library(tidyr)
library(knitr)
library(ggplot2)
library(Rcpp)
library(cowplot) # 新增：用于提取和保存图例

# --- 全局设置 ---
i <- 1 # 你的模拟设置编号 对应的case
n_values <- c(1000) # 要循环的n值
q_track_alpha <- c(1, 2)
q_track_beta <- c(1, 2)
q_track_gamma <- c(1)
K <- 100
t_grid <- round(seq(1/K, 1, 1/K), 2)
coverage_points <- c(0.3, 0.5, 0.7)

# --- 编译 C++ 代码 (只需一次) ---
sourceCpp(file = paste0("src/SimSet", i, ".cpp"))

# ===================================================================
# 1. HELPER FUNCTIONS - 提取指标、绘图、生成表格
# ===================================================================

#' @title Extract Metrics from Simulation Results
extract_metrics <- function(results_hat, results_var, true_func, t_grid, coverage_points) {
  hat_matrix <- do.call(rbind, results_hat)
  var_matrix <- do.call(rbind, results_var)
  mise <- mean(((hat_matrix - matrix(true_func, nrow = nrow(hat_matrix), ncol = length(t_grid), byrow = TRUE))^2), na.rm = TRUE)
  coverage_indices <- which(t_grid %in% coverage_points)
  coverage_rates <- numeric(length(coverage_indices))
  ci_lengths <- numeric(length(coverage_indices))
  if (length(coverage_indices) > 0) {
    for (j in 1:length(coverage_indices)) {
      idx <- coverage_indices[j]
      # 这里的CI长度和覆盖率是基于估计的方差计算的
      coverage_rates[j] <- mean(abs(hat_matrix[, idx] - true_func[idx]) / sqrt(var_matrix[, idx]) < 1.96, na.rm = TRUE)
      ci_lengths[j] <- mean(2 * 1.96 * sqrt(var_matrix[, idx]), na.rm = TRUE)
    }
  }
  return(list(MISE = mise, coverage_rates = coverage_rates, ci_lengths = ci_lengths, coverage_points = t_grid[coverage_indices]))
}


#' @title Generate and Save a Plot
#' @description This function generates a plot showing the mean estimate, true function,
#'              empirical 95% CI bands, and estimated 95% CI bands.
generate_and_save_plot <- function(hat_list, var_list, true_func, t_grid, ylab_expr, analysis_type, file_path) {
  hat_matrix <- do.call(rbind, hat_list)
  var_matrix <- do.call(rbind, var_list) # 需要 var_matrix 来计算估计的CI
  
  plot_df <- data.frame(
    t = t_grid,
    mean_est = colMeans(hat_matrix, na.rm = TRUE),
    # 经验的95% CI (基于估计值的2.5%和97.5%分位数)
    lower_ci_empirical = apply(hat_matrix, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper_ci_empirical = apply(hat_matrix, 2, quantile, probs = 0.975, na.rm = TRUE),
    # 真实函数值
    true_val = true_func
  )
  
  # 计算估计的95% CI (平均估计值 +/- 1.96 * sqrt(平均估计方差))
  mean_var_est <- colMeans(var_matrix, na.rm = TRUE)
  
  # 确保平均方差非负，避免sqrt(负数)产生NaN
  mean_var_est[mean_var_est < 0] <- 0 
  
  plot_df <- plot_df %>%
    mutate(
      lower_ci_estimated = mean_est - 1.96 * sqrt(mean_var_est),
      upper_ci_estimated = mean_est + 1.96 * sqrt(mean_var_est)
    )
  
  # 定义线条和填充的颜色 (美观且易描述的颜色)
  line_colors <- c("Estimate" = "black", "True" = "red")
  fill_colors <- c("Empirical 95% CI" = "blue", "Estimated 95% CI" = "green") 
  
  p <- ggplot(plot_df, aes(x = t)) +
    # 绘制估计的95% CI bands (稍低的透明度，使其在经验CI下方或与经验CI区分)
    geom_ribbon(aes(ymin = lower_ci_estimated, ymax = upper_ci_estimated, fill = "Estimated 95% CI"), alpha = 0.25) +
    # 绘制经验的95% CI bands (稍高的透明度，使其在估计CI上方)
    geom_ribbon(aes(ymin = lower_ci_empirical, ymax = upper_ci_empirical, fill = "Empirical 95% CI"), alpha = 0.25) +
    # 绘制平均估计线
    geom_line(aes(y = mean_est, color = "Estimate")) +
    # 绘制真实函数线
    geom_line(aes(y = true_val, color = "True"), linetype = "dashed") +
    
    # 手动设置线条颜色和填充颜色，并移除图例标题
    scale_color_manual(name = NULL, values = line_colors) +
    scale_fill_manual(name = NULL, values = fill_colors) +
    
    labs(x = "Time (t)", y = ylab_expr) +
    theme_bw() +
    theme(
      legend.position = "none" # 移除所有图例
    )
  
  if (analysis_type == "smooth") {
    p <- p + coord_cartesian(xlim = c(0.1, 0.9))
  }
  ggsave(filename = file_path, plot = p, width = 3, height = 2.5, device = "pdf")
}


#' @title Create and Print Summary Tables (MODIFIED for RMISE)
create_summary_tables <- function(data, table_title_prefix) {
  if (nrow(data) == 0) {
    cat(paste("\nNo data available to generate tables for:", table_title_prefix, "\n"))
    return()
  }
  
  # --- Table 1: Coverage Rate (CI Length) ---
  median_ci <- median(data$ci_length, na.rm = TRUE)
  ci_multiplier <- 1
  caption_note <- ""
  
  if (median_ci > 0 && median_ci < 0.1) {
    exponent <- floor(log10(median_ci))
    if (exponent != 0) {
      ci_multiplier <- 10^(-exponent)
      caption_note <- paste0(" (CI Lengths scaled by 10^", -exponent, ")")
    }
  }
  
  table1_data <- data %>%
    mutate(
      scaled_ci = ci_length * ci_multiplier,
      formatted_result = sprintf("%.1f (%.2f)", coverage_rate * 100, scaled_ci)
    ) %>%
    dplyr::select(n, parameter, t, formatted_result) %>%
    pivot_wider(names_from = t, values_from = formatted_result, names_prefix = "t=") %>%
    arrange(n, parameter)
  
  cat(paste0("\n\n--- ", table_title_prefix, ": Table 1: Coverage Rate (%) and Mean CI Length ---\n"))
  print(kable(table1_data, 
              caption = paste0(table_title_prefix, " - Coverage Rate (%) and Mean CI Length", caption_note),
              align = 'c'))
  
  # --- Table 2: RMISE (MODIFIED from MISE) ---
  table2_data <- data %>%
    dplyr::select(n, parameter, MISE) %>%
    distinct() %>%
    mutate(RMISE = sqrt(MISE)) %>% # Calculate Root MISE
    dplyr::select(-MISE) %>% # Remove the original MISE column
    pivot_wider(names_from = parameter, values_from = RMISE) %>%
    arrange(n)
  
  cat(paste0("\n\n--- ", table_title_prefix, ": Table 2: Root Mean Integrated Squared Error (RMISE) ---\n"))
  print(kable(table2_data, 
              caption = paste(table_title_prefix, "- Root MISE (RMISE)"),
              digits = 5, # Adjust digits for RMISE scale
              align = 'c'))
}


# ===================================================================
# 2. MAIN ANALYSIS LOOP
# ===================================================================

all_results_list <- list()
cat("Starting analysis for all n values and analysis types...\n")

for (n_val in n_values) {
  cat(paste("\nProcessing n =", n_val, "...\n"))
  load_path <- paste0("result/test_simu_", i, "/n=", n_val, ".rdata")
  if (!file.exists(load_path)) {
    warning(paste("File not found, skipping:", load_path))
    next
  }
  load(file = load_path)
  
  for (analysis_type in c("smooth", "integrated")) {
    cat(paste("  Analyzing", analysis_type, "functions...\n"))
    
    if (analysis_type == "smooth") {
      param_config <- list(
        list(name = "alpha", type = "alpha", indices = q_track_alpha, true_func = function(t, q) unlist(lapply(t, fs, i = q, n = n_val))),
        list(name = "beta", type = "beta", indices = q_track_beta, true_func = function(t, q) unlist(lapply(t, fr, j = q, n = n_val))),
        list(name = "gamma", type = "gamma", indices = q_track_gamma, true_func = function(t, q) unlist(lapply(t, fg, n = n_val)))
      )
    } else { # integrated
      param_config <- list(
        list(name = "Alpha", type = "alpha", indices = q_track_alpha, true_func = function(t, q) fsI_cpp(t, i = q, n = n_val)),
        list(name = "Beta", type = "beta", indices = q_track_beta, true_func = function(t, q) frI_cpp(t, j = q, n = n_val)),
        list(name = "Gamma", type = "gamma", indices = q_track_gamma, true_func = function(t, q) fgI_cpp(t, n = n_val))
      )
    }
    
    for (p in param_config) {
      for (q in p$indices) {
        param_name_full <- paste0(p$name, "_", q)
        q_name <- paste0("q_", q)
        
        hat_list <- results[[p$type]][[analysis_type]]$hat[[q_name]]
        var_list <- results[[p$type]][[analysis_type]]$var[[q_name]]
        true_values <- p$true_func(t_grid, q)
        
        metrics <- extract_metrics(hat_list, var_list, true_values, t_grid, coverage_points)
        result_df <- tibble(
          analysis_type = analysis_type, n = n_val, parameter = param_name_full,
          MISE = metrics$MISE, t = metrics$coverage_points,
          coverage_rate = metrics$coverage_rates, ci_length = metrics$ci_lengths
        )
        all_results_list[[length(all_results_list) + 1]] <- result_df
        
        plot_ylab_expr <- switch(p$name,
                                 "alpha" = bquote(hat(alpha)[.(q)](t)),
                                 "beta"  = bquote(hat(beta)[.(q)](t)),
                                 "gamma" = bquote(hat(gamma)[.(q)](t)),
                                 "Alpha" = bquote(hat(A)[.(q)](t)),
                                 "Beta"  = bquote(hat(B)[.(q)](t)),
                                 "Gamma" = bquote(hat(Gamma)[.(q)](t))
        )
        
        plot_filepath <- paste0("result/test_simu_", i, "/", param_name_full, "_n", n_val, ".pdf")
        
        generate_and_save_plot(
          hat_list = hat_list,
          var_list = var_list,
          true_func = true_values,
          t_grid = t_grid,
          ylab_expr = plot_ylab_expr,
          analysis_type = analysis_type,
          file_path = plot_filepath
        )
      }
    }
  }
}
cat("\nAnalysis and plotting complete. Generating summary tables...\n")

# ===================================================================
# 3. FINAL TABLE GENERATION
# ===================================================================
results_df <- bind_rows(all_results_list)

create_summary_tables(
  data = results_df %>% filter(analysis_type == "smooth"),
  table_title_prefix = "Smooth Functions (alpha, beta, gamma)"
)

create_summary_tables(
  data = results_df %>% filter(analysis_type == "integrated"),
  table_title_prefix = "Integrated Functions (Alpha, Beta, Gamma)"
)
