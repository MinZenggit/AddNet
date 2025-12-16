# 文件名: nonParametric_diff_flexible3.R (已更新)
# The C++ files (compute_flexible2.cpp, create_tensors_flexible.cpp)
# should be in the working directory and will be sourced automatically by Rcpp
# if you use sourceCpp, or handled by your project's compilation setup.

#' Flexible Non-parametric Difference Estimator with Time-Varying Covariates
#'
#' Computes non-parametric estimates using truncation and (optionally) post-smoothing.
#' This version supports covariates Z that change at specified time points.
#' The smoothing is performed directly on the truncated point estimates.
#'
#' @param data A data frame with event data (sender, receiver, time).
#' @param z_variants A 4D array (n x n x p x Kz) of Kz distinct covariate sets.
#'                   Set to NULL if p = 0.
#' @param z_times A numeric vector of length Kz specifying the normalized start times [0, 1) 
#'                for each covariate set in z_variants. Must be sorted.
#' @param n The number of nodes.
#' @param p The number of covariates.
#' @param K The number of time points/bins.
#' @param smooth A logical value. If `TRUE`, post-smoothing results are computed.
#' @param h1 The bandwidth for smoothing alpha and beta estimates. Required if `smooth = TRUE`.
#' @param h2 The bandwidth for smoothing gamma estimates. Required if `smooth = TRUE`.
#'
#' @return A list containing results for `truncate` and optionally `smooth` methods.
#'
nonParametric_diff_flexible <- function(data, z_variants, z_times, n, p, K, 
                                        smooth = FALSE, 
                                        h1 = NULL, h2 = NULL) {
  
  if (smooth && (is.null(h1) || is.null(h2))) {
    stop("`h1` and `h2` must be provided when `smooth = TRUE`.")
  }
  
  # --- 1. Input Validation for Time-Varying Z ---
  if (p > 0) {
    if (is.null(z_variants) || is.null(z_times)) {
      stop("z_variants and z_times must be provided when p > 0.")
    }
    Kz <- dim(z_variants)[4]
    if (length(z_times) != Kz) {
      stop("Length of z_times must match the 4th dimension of z_variants (Kz).")
    }
  } else {
    Kz <- 1 # A single "empty" Z set
  }
  
  # --- 2. Time Normalization and Event Tensor Creation ---
  # We only ever need the truncated event counts from C++.
  # The `smooth` flag now controls post-processing in R.
  t_start <- data[1, ncol(data)]
  t_end <- data[nrow(data), ncol(data)]
  t_total <- t_end - t_start
  data$t_norm <- (data[, ncol(data)] - t_start) / t_total
  
  time_points <- seq(from = 1/(2*K), to = 1 - 1/(2*K), by = 1/K)
  
  event_tensors <- create_event_tensors_flexible_cpp(
    senders = data[, 1],
    receivers = data[, 2],
    event_times = data$t_norm,
    n = n, K = K,
    bandwidth = 0.1, # Dummy value, not used for truncate
    create_kernel = FALSE # IMPORTANT: We never create the kernel tensor in C++ anymore
  )
  
  # --- 3. PRE-COMPUTATION OF Z-DERIVATIVES (UNCHANGED) ---
  precomputed_z_list <- vector("list", Kz)
  if (p > 0) {
    for (j in 1:Kz) {
      z_current <- z_variants[,,,j, drop = FALSE]
      diff2_z_mat <- (z_current[2:n, 2:n, , ,drop = FALSE] - z_current[1:(n-1), 2:n, , ,drop = FALSE] - 
                        z_current[2:n, 1:(n-1), , ,drop = FALSE] + z_current[1:(n-1), 1:(n-1), , ,drop = FALSE])
      dim(diff2_z_mat) <- c((n-1)^2, p)
      XtX <- t(diff2_z_mat) %*% diff2_z_mat
      if (rcond(XtX) < .Machine$double.eps) {
        stop(paste("Covariate matrix for Z-variant", j, "is singular."))
      }
      inv_XtX <- solve(XtX)
      diff2_z_core <- inv_XtX %*% t(diff2_z_mat)
      diff2_z_cube <- array(diff2_z_mat, dim = c(n-1, n-1, p))
      z_diff_alpha <- z_current[2:n, , , ,drop=FALSE] - z_current[1:(n-1), , , ,drop=FALSE]
      z_diff_beta <- z_current[, 2:n, , ,drop=FALSE] - z_current[, 1:(n-1), , ,drop=FALSE]
      precomputed_z_list[[j]] <- list(
        z_k = array(z_current, dim=c(n,n,p)), inv_XtX = inv_XtX, diff2_z_core = diff2_z_core,
        diff2_z_cube = diff2_z_cube, z_diff_alpha = array(z_diff_alpha, dim=c(n-1,n,p)),
        z_diff_beta = array(z_diff_beta, dim=c(n,n-1,p))
      )
    }
  } else {
    precomputed_z_list[[1]] <- list(
      z_k = array(0, dim=c(n,n,0)), inv_XtX = matrix(0,0,0), diff2_z_core = matrix(0,0,(n-1)^2),
      diff2_z_cube = array(0, dim=c(n-1,n-1,0)), z_diff_alpha = array(0, dim=c(n-1,n,0)),
      z_diff_beta = array(0, dim=c(n,n-1,0))
    )
  }
  
  # --- 4. Create Mapping from Time Slice (k) to Z-set (j) (UNCHANGED) ---
  if (p > 0) {
    z_map_indices <- findInterval(time_points, z_times, rightmost.closed = TRUE)
    z_map_indices[z_map_indices == 0] <- 1
  } else {
    z_map_indices <- rep(1, K)
  }
  
  # --- 5. Call C++ Computation Engine (Now only for truncate results) ---
  raw_results <- compute_results_flexible_cpp(
    n = n, p = p,
    event_counts_truncate = event_tensors$truncate,
    event_counts_kernel_nullable = NULL, # We pass NULL for the kernel counts
    precomputed_z_list = precomputed_z_list,
    z_map_indices = z_map_indices
  )
  
  # --- 6. Post-processing and Result Assembly ---
  
  # Helper for integration (UNCHANGED)
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
  
  # Assemble truncate results (UNCHANGED)
  results_truncate <- list(
    pointed = raw_results$truncate,
    integrated = integrate_results(raw_results$truncate, K, p)
  )
  
  final_output <- list(
    time_points = time_points,
    truncate = results_truncate
  )
  
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
    
    # Add confidence intervals to the pointed results
    pointed_smooth$alpha_CI <- array(c(pointed_smooth$alpha_hat - 1.96 * sqrt(pointed_smooth$alpha_var),
                                       pointed_smooth$alpha_hat + 1.96 * sqrt(pointed_smooth$alpha_var)),
                                     dim = c(n, K, 2))
    pointed_smooth$beta_CI <- array(c(pointed_smooth$beta_hat - 1.96 * sqrt(pointed_smooth$beta_var),
                                      pointed_smooth$beta_hat + 1.96 * sqrt(pointed_smooth$beta_var)),
                                    dim = c(n, K, 2))
    if (p > 0) {
      gamma_se <- t(apply(pointed_smooth$gamma_var, 3, function(mat) sqrt(diag(mat))))
      pointed_smooth$gamma_CI <- array(c(pointed_smooth$gamma_hat - 1.96 * gamma_se,
                                         pointed_smooth$gamma_hat + 1.96 * gamma_se),
                                       dim = c(p, K, 2))
    }
    
    final_output$smooth <- list(
      pointed = pointed_smooth,
      integrated = integrate_results(pointed_smooth, K, p)
    )
  }
  
  return(final_output)
}
