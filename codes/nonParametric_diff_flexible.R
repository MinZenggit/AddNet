

#' Flexible Non-parametric Difference Estimator
#'
#' Computes non-parametric estimates using the truncation method by default.
#' If `smooth = TRUE`, it also computes estimates using the kernel smoothing method.
#'
#' @param data A data frame with event data (sender, receiver, time).
#' @param zij A 4D array of covariates (n x n x p x K).
#' @param n The number of nodes.
#' @param p The number of covariates.
#' @param K The number of time points/bins.
#' @param smooth A logical value. If `TRUE`, kernel smoothing results are also computed. Defaults to `FALSE`.
#' @param bandwidth The bandwidth `h` for kernel smoothing. Required if `smooth = TRUE`.
#'
#' @return A list containing results. The `truncate` element is always present.
#'   The `kernel` element is only present if `smooth = TRUE`.
#'
nonParametric_diff_flexible <- function(data, zij, n, p, K, 
                                        smooth = FALSE, 
                                        bandwidth = NULL) {
  
  if (smooth && is.null(bandwidth)) {
    stop("`bandwidth` must be provided when `smooth = TRUE`.")
  }
  
  # --- 1. 通用预处理 ---
  t_start <- data[1, ncol(data)]
  t_end <- data[nrow(data), ncol(data)]
  t_total <- t_end - t_start
  data$t_norm <- (data[, ncol(data)] - t_start) / t_total
  
  time_points <- seq(from = 1/(2*K), to = 1 - 1/(2*K), by = 1/K)
  
  # if (p > 0 && abs(sum(zij[,,,1])) > 1e-9) {
  #   warning("Sum of covariates `zij` is not zero. Variance formula for lambda0 might be inaccurate.")
  # }
  
  if (p > 0) {
    z_constant <- array(zij[, , , 1], dim = c(n, n, p))
    diff2_z_mat <- (z_constant[2:n, 2:n, , drop = FALSE] - z_constant[1:(n-1), 2:n, , drop = FALSE] - 
                      z_constant[2:n, 1:(n-1), , drop = FALSE] + z_constant[1:(n-1), 1:(n-1), , drop = FALSE])
    dim(diff2_z_mat) <- c((n-1)^2, p)
    XtX <- t(diff2_z_mat) %*% diff2_z_mat
    if (rcond(XtX) < .Machine$double.eps) stop("Covariate matrix is singular.")
    inv_XtX <- solve(XtX)
    diff2_z_core <- inv_XtX %*% t(diff2_z_mat)
    diff2_z_cube <- array(diff2_z_mat, dim = c(n-1, n-1, p))
  } else {
    z_constant <- array(0, dim = c(n, n, 0))
    diff2_z_core <- matrix(0, nrow = 0, ncol = (n-1)^2)
    inv_XtX <- matrix(0, nrow = 0, ncol = 0)
    diff2_z_cube <- array(0, dim = c(n-1, n-1, 0))
  }
  
  # --- 2. 调用C++生成事件张量 (按需) ---
  # message(paste("Step 1: Creating event tensor(s)...", if(smooth) "(Truncate + Kernel)" else "(Truncate only)"))
  event_tensors <- create_event_tensors_flexible_cpp(
    senders = data[, 1],
    receivers = data[, 2],
    event_times = data$t_norm,
    n = n, K = K,
    bandwidth = if(is.null(bandwidth)) 0.1 else bandwidth, # Pass a dummy value if not used
    create_kernel = smooth
  )
  
  # --- 3. 调用C++计算参数 (按需) ---
  # message(paste("Step 2: Computing parameters...", if(smooth) "(Truncate + Kernel)" else "(Truncate only)"))
  raw_results <- compute_results_flexible_cpp(
    n = n, p = p,
    event_counts_truncate = event_tensors$truncate,
    event_counts_kernel_nullable = event_tensors$kernel, # This can be NULL
    z_constant = z_constant,
    diff2_z_core = diff2_z_core,
    inv_XtX = inv_XtX,
    diff2_z_cube = diff2_z_cube
  )
  
  # --- 4. 后处理和结果组装 ---
  # message("Step 3: Post-processing and assembling results...")
  
  # Helper function for integration
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
      integrated$gamma_var <- array(0, dim = c(p, p, K))
      integrated$gamma_var[,,1] <- pointed_results$gamma_var[,,1]
      for (k in 2:K) {
        integrated$gamma_var[,,k] <- integrated$gamma_var[,,k-1] + pointed_results$gamma_var[,,k]
      }
      integrated$gamma_var <- integrated$gamma_var * (integration_step^2)
    } else {
      integrated$gamma_hat <- NULL
      integrated$gamma_var <- NULL
    }
    return(integrated)
  }
  
  # Process Truncate results (always)
  results_truncate <- list(
    pointed = raw_results$truncate,
    integrated = integrate_results(raw_results$truncate, K, p)
  )
  
  # Initialize final output list
  final_output <- list(
    time_points = time_points,
    truncate = results_truncate
  )
  
  # Process Kernel results (conditionally)
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

