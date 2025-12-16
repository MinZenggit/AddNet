# nonParametric_diff_flexible_tv.R

#' Flexible Non-parametric Difference Estimator for Time-Varying Covariates
#'
#' Computes non-parametric estimates using a sparse backend, supporting covariates
#' that are piecewise constant over time.
#'
#' @param data A data frame with event data (sender, receiver, time).
#' @param z_list A list for time-varying covariates. Required if p > 0.
#'   It must contain:
#'   - `$times`: A numeric vector of time points where covariates change.
#'              These should be on the same scale as event times in `data`.
#'              For K intervals, K-1 time points are needed. e.g., c(0.1, 0.3) for 3 intervals.
#'   - `$covariates`: A list of sparse data frames. Each data frame corresponds
#'                   to a time interval and must contain `sender`, `receiver`,
#'                   and covariate columns.
#' @param n The number of nodes.
#' @param p The number of covariates.
#' @param K The number of time points/bins for estimation.
#' @param smooth A logical value. If `TRUE`, post-smoothing results are computed.
#' @param h1 Bandwidth for smoothing alpha/beta.
#' @param h2 Bandwidth for smoothing gamma.

#' @return A list containing truncated estimation results.
#'
nonParametric_diff_flexible_tv <- function(data, z_list = NULL, n, p, K, 
                                           smooth = FALSE, 
                                           h1 = NULL, h2 = NULL) {
  
  # --- 0. 输入验证 (已更新) ---
  if (p > 0) {
    if (is.null(z_list)) {
      stop("`z_list` must be provided when `p > 0`.")
    }
    if (!is.list(z_list) || !all(c("times", "covariates") %in% names(z_list))) {
      stop("`z_list` must be a list containing 'times' and 'covariates'.")
    }
    if (length(z_list$times) != length(z_list$covariates) - 1) {
      stop("Length of `z_list$times` must be one less than the length of `z_list$covariates`.")
    }
    if (!is.list(z_list$covariates)) {
      stop("`z_list$covariates` must be a list of data frames.")
    }
  }
  if (p == 0 && !is.null(z_list)) {
    warning("`p` is 0 but `z_list` is provided. It will be ignored.")
    z_list = NULL
  }
  
  # --- 1. 通用预处理 (无变化) ---
  t_start <- min(data$t)
  t_end <- max(data$t)
  t_total <- t_end - t_start
  # data$t_norm <- (data$t - t_start) / t_total
  
  time_points <- seq(from = 1/(2*K), to = 1 - 1/(2*K), by = 1/K)
  
  # --- 2. 协变量格式化 (已更新) ---
  z_list_cpp <- NULL
  z_times_norm <- NULL
  
  if (p > 0) {
    # 将每个协变量数据帧转换为C++所需的格式 (带有一个矩阵列)
    z_list_cpp <- lapply(z_list$covariates, function(z_df) {
      z_values_mat <- as.matrix(z_df[, -(1:2), drop = FALSE])
      z_df_cpp <- data.frame(
        sender = z_df$sender,
        receiver = z_df$receiver
      )
      z_df_cpp$values <- z_values_mat
      return(z_df_cpp)
    })
    
    # 对协变量变化的时间点进行归一化
    z_times_norm <- (z_list$times - t_start) / t_total
  }
  
  # --- 3. 调用C++生成稀疏事件列表 (无变化) ---
  # message("Step 1/3: Creating sparse event lists...")
  # 注意：这里我们只创建 truncate 部分，因为 smooth 部分已移除
  sparse_event_dfs <- create_event_tensors_sparse_cpp(
    senders = data$s,
    receivers = data$r,
    event_times = data$t,
    n = n, K = K,
    bandwidth = 0.1, # 占位符，C++内部不使用
    create_kernel = FALSE # 明确告知不创建kernel数据
  )
  
  # --- 4. 调用C++进行核心计算 (已更新为新的C++函数) ---
  # message("Step 2/3: Computing parameters from sparse data...")
  # 使用新的C++函数，它能处理随时间变化的协变量
  raw_results <- compute_results_sparse_cpp_tv(
    n = n, p = p, K = K,
    events_df = sparse_event_dfs$truncate,
    z_list_nullable = z_list_cpp,
    z_times_nullable = z_times_norm
  )
  
  # --- 5. 后处理和结果组装 (已简化，移除smooth部分) ---
  # message("Step 3/3: Post-processing and assembling results...")
  
  integrate_results <- function(pointed_results, K, p) {
    integration_step <- 1.0 / K
    integrated <- list()
    
    integrated$lambda0_hat <- cumsum(pointed_results$lambda0_hat) * integration_step
    integrated$lambda0_var <- cumsum(pointed_results$lambda0_var) * (integration_step^2)
    integrated$alpha_hat <- t(apply(pointed_results$alpha_hat, 1, cumsum)) * integration_step
    integrated$alpha_var <- t(apply(pointed_results$alpha_var, 1, cumsum)) * (integration_step^2)
    integrated$beta_hat <- t(apply(pointed_results$beta_hat, 1, cumsum)) * integration_step
    integrated$beta_var <- t(apply(pointed_results$beta_var, 1, cumsum)) * (integration_step^2)
    
    if (p > 0) {
      integrated$gamma_hat <- t(apply(pointed_results$gamma_hat, 1, cumsum)) * integration_step
      
      cumm_gamma_var <- apply(pointed_results$gamma_var, c(1, 2), cumsum)
      if (K == 1) {
        dim(cumm_gamma_var) <- c(1, p, p)
      }
      cumm_gamma_var_permuted <- aperm(cumm_gamma_var, c(2, 3, 1))
      integrated$gamma_var <- cumm_gamma_var_permuted * (integration_step^2)
      
    } else {
      integrated$gamma_hat <- NULL
      integrated$gamma_var <- NULL
    }
    return(integrated)
  }
  
  # 注意：raw_results现在直接就是 truncate 的结果
  results_truncate <- list(
    pointed = raw_results,
    integrated = integrate_results(raw_results, K, p)
  )
  
  final_output <- list(
    time_points = time_points,
    truncate = results_truncate
  )
  
  # message("Done.")
  # --- 7. NEW: POST-SMOOTHING BLOCK ---
  if (smooth) {
    # Get raw estimates and variances from the truncate results
    raw_est <- final_output$truncate$pointed
    raw_var <- final_output$truncate$pointed # Variances are stored with point estimates in C++ output
    
    # Define kernel function
    ker_epan <- function(u) { 0.75 * (1 - u^2) * (abs(u) <= 1) }
    
    # Initialize storage for smoothed results
    alpha_hat_smooth <- matrix(0, nrow = n, ncol = K)
    beta_hat_smooth  <- matrix(0, nrow = n, ncol = K)
    V_alpha_smooth   <- matrix(0, nrow = n, ncol = K)
    V_beta_smooth    <- matrix(0, nrow = n, ncol = K)
    
    if (p > 0) {
      gamma_hat_smooth <- matrix(0, nrow = p, ncol = K)
      V_gamma_smooth   <- array(0, dim = c(p, p, K))
    }
    
    # Loop over each time point to calculate smoothed value
    for (k in 1:K) {
      t_k <- time_points[k]
      
      # --- Smoothing for alpha and beta (bandwidth h1) ---
      weights_h1 <- ker_epan((time_points - t_k) / h1)
      sum_weights_h1 <- sum(weights_h1)
      
      if (sum_weights_h1 > 1e-9) {
        W_k_h1 <- weights_h1 / sum_weights_h1
        
        # Smooth point estimates
        alpha_hat_smooth[, k] <- raw_est$alpha_hat %*% W_k_h1
        beta_hat_smooth[, k]  <- raw_est$beta_hat %*% W_k_h1
        
        # Smooth variance using Var(Sy) = S Var(y) S' approx sum(w_i^2 * Var(y_i))
        V_alpha_smooth[, k] <- raw_var$alpha_var %*% (W_k_h1^2)
        V_beta_smooth[, k]  <- raw_var$beta_var %*% (W_k_h1^2)
      } else { # Fallback if bandwidth is too small
        alpha_hat_smooth[, k] <- raw_est$alpha_hat[, k]
        beta_hat_smooth[, k]  <- raw_est$beta_hat[, k]
        V_alpha_smooth[, k]   <- raw_var$alpha_var[, k]
        V_beta_smooth[, k]    <- raw_var$beta_var[, k]
      }
      
      # --- Smoothing for gamma (bandwidth h2) ---
      if (p > 0) {
        weights_h2 <- ker_epan((time_points - t_k) / h2)
        sum_weights_h2 <- sum(weights_h2)
        
        if (sum_weights_h2 > 1e-9) {
          W_k_h2 <- weights_h2 / sum_weights_h2
          
          # Smooth point estimates
          gamma_hat_smooth[, k] <- raw_est$gamma_hat %*% W_k_h2
          
          # Smooth variance
          # For each (l,m) entry of the covariance matrix, smooth it
          for (l in 1:p) {
            for (m in 1:p) {
              raw_gamma_cov_lm <- raw_var$gamma_var[l, m, ]
              V_gamma_smooth[l, m, k] <- sum(raw_gamma_cov_lm * (W_k_h2^2))
            }
          }
        } else { # Fallback
          gamma_hat_smooth[, k] <- raw_est$gamma_hat[, k]
          V_gamma_smooth[,, k]  <- raw_var$gamma_var[,, k]
        }
      }
    }
    
    # --- Assemble smoothed results ---
    pointed_smooth <- list(
      alpha_hat = alpha_hat_smooth,
      beta_hat = beta_hat_smooth,
      alpha_var = V_alpha_smooth,
      beta_var = V_beta_smooth,
      # Dummy lambda0 for now, can be added if needed
      lambda0_hat = raw_est$lambda0_hat, 
      lambda0_var = raw_var$lambda0_var
    )
    
    if (p > 0) {
      pointed_smooth$gamma_hat <- gamma_hat_smooth
      pointed_smooth$gamma_var <- V_gamma_smooth
    }
    final_output$smooth <- list(
      pointed = pointed_smooth,
      integrated = integrate_results(pointed_smooth, K, p)
    )
  }
  
  return(final_output)
}
