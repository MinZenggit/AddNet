# 确保已安装和加载必要的包
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}
library(Matrix)

#' Flexible Non-parametric Difference Estimator (High-Performance Sparse Version)
#'
#' Computes non-parametric estimates using a fully sparse backend to handle very large n.
#'
#' @param data A data frame with event data (sender, receiver, time).
#' @param z_df A data frame for time-invariant covariates in sparse "long" format. 
#'   Must contain columns `sender`, `receiver`, and one column for each covariate.
#'   Only non-zero covariate pairs need to be included. Can be NULL if p=0.
#' @param n The number of nodes.
#' @param p The number of covariates.
#' @param K The number of time points/bins.
#' @param smooth A logical value. If `TRUE`, kernel smoothing results are also computed.
#' @param bandwidth The bandwidth `h` for kernel smoothing. Required if `smooth = TRUE`.
#'
#' @return A list containing results.
#'
nonParametric_diff_flexible <- function(data, z_df = NULL, n, p, K, 
                                        smooth = FALSE, 
                                        bandwidth = NULL) {
  
  # --- 0. 输入验证 ---
  if (smooth && is.null(bandwidth)) {
    stop("`bandwidth` must be provided when `smooth = TRUE`.")
  }
  if (p > 0 && is.null(z_df)) {
    stop("`z_df` must be provided when `p > 0`.")
  }
  if (p == 0 && !is.null(z_df)) {
    warning("`p` is 0 but `z_df` is provided. It will be ignored.")
    z_df = NULL
  }
  
  # --- 1. 通用预处理 ---
  t_start <- data[1, ncol(data)]
  t_end <- data[nrow(data), ncol(data)]
  t_total <- t_end - t_start
  data$t_norm <- (data[, ncol(data)] - t_start) / t_total
  
  time_points <- seq(from = 1/(2*K), to = 1 - 1/(2*K), by = 1/K)
  
  # --- 2. 协变量格式化 ---
  # 将协变量数据帧转换为C++所需的格式
  if (p > 0) {
    # 将所有协变量列合并到一个名为 "values" 的矩阵列中
    z_values_mat <- as.matrix(z_df[, -(1:2), drop = FALSE])
    z_df_cpp <- data.frame(
      sender = z_df$sender,
      receiver = z_df$receiver
    )
    z_df_cpp$values <- z_values_mat
  } else {
    z_df_cpp <- NULL
  }
  
  # --- 3. 调用C++生成稀疏事件列表 ---
  # message("Step 1/3: Creating sparse event lists...")
  sparse_event_dfs <- create_event_tensors_sparse_cpp(
    senders = data[, 1],
    receivers = data[, 2],
    event_times = data$t_norm,
    n = n, K = K,
    bandwidth = if(is.null(bandwidth)) 0.1 else bandwidth,
    create_kernel = smooth
  )
  
  # --- 4. 调用C++进行核心计算 (使用新的稀疏接口) ---
  # message("Step 2/3: Computing parameters from sparse data (this may take a while)...")
  raw_results <- compute_results_sparse_cpp(
    n = n, p = p, K = K,
    events_df = sparse_event_dfs$truncate,
    z_df_nullable = z_df_cpp,
    compute_kernel = smooth,
    events_kernel_df_nullable = sparse_event_dfs$kernel
  )
  
  # --- 5. 后处理和结果组装 (逻辑不变) ---
  # message("Step 3/3: Post-processing and assembling results...")
  
  integrate_results <- function(pointed_results, K, p) {
    integration_step <- 1.0 / K
    integrated <- list()
    
    # lambda0, alpha, beta 的计算保持不变
    integrated$lambda0_hat <- cumsum(pointed_results$lambda0_hat) * integration_step
    integrated$lambda0_var <- cumsum(pointed_results$lambda0_var) * (integration_step^2)
    integrated$alpha_hat <- t(apply(pointed_results$alpha_hat, 1, cumsum)) * integration_step
    integrated$alpha_var <- t(apply(pointed_results$alpha_var, 1, cumsum)) * (integration_step^2)
    integrated$beta_hat <- t(apply(pointed_results$beta_hat, 1, cumsum)) * integration_step
    integrated$beta_var <- t(apply(pointed_results$beta_var, 1, cumsum)) * (integration_step^2)
    
    if (p > 0) {
      # gamma_hat 的计算保持不变
      integrated$gamma_hat <- t(apply(pointed_results$gamma_hat, 1, cumsum)) * integration_step
      
      # ************************** BUG FIX STARTS HERE **************************
      # 废弃手写的 for 循环，改用与 alpha/beta_var 一致的 apply 方法
      
      # 1. 对 p x p x K 的数组，沿着第3个维度(K)进行累积求和
      #    apply 的结果是一个 K x p x p 的数组
      cumm_gamma_var <- apply(pointed_results$gamma_var, c(1, 2), cumsum)
      
      # 2. 如果 K=1, apply 会降维成 p x p 矩阵, 需要手动加回第三维
      if (K == 1) {
        dim(cumm_gamma_var) <- c(1, p, p)
      }
      
      # 3. 使用 aperm 将维度从 K x p x p 转回 p x p x K
      cumm_gamma_var_permuted <- aperm(cumm_gamma_var, c(2, 3, 1))
      
      # 4. 乘以积分步长的平方
      integrated$gamma_var <- cumm_gamma_var_permuted * (integration_step^2)
      # ************************** BUG FIX ENDS HERE ****************************
      
    } else {
      integrated$gamma_hat <- NULL
      integrated$gamma_var <- NULL
    }
    return(integrated)
  }
  
  
  results_truncate <- list(
    pointed = raw_results$truncate,
    integrated = integrate_results(raw_results$truncate, K, p)
  )
  
  final_output <- list(
    time_points = time_points,
    truncate = results_truncate
  )
  
  if (smooth) {
    kernel_var_scaler <- 0.28209479177 / bandwidth
    pointed_kernel <- raw_results$kernel
    pointed_kernel$lambda0_var <- pointed_kernel$lambda0_var * kernel_var_scaler
    pointed_kernel$alpha_var   <- pointed_kernel$alpha_var * kernel_var_scaler
    pointed_kernel$beta_var    <- pointed_kernel$beta_var * kernel_var_scaler
    if (p > 0) {
      pointed_kernel$gamma_var <- pointed_kernel$gamma_var * kernel_var_scaler
    }
    
    final_output$kernel <- list(
      pointed = pointed_kernel,
      integrated = integrate_results(pointed_kernel, K, p)
    )
  }
  
  # message("Done.")
  return(final_output)
}
