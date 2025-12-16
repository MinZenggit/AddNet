# 确保 Matrix 包已加载
if (!require("Matrix")) {
  install.packages("Matrix")
}
library(Matrix)

#' @title Calculate Dynamic Pair-wise Network Covariates (v3)
#' @description Efficiently calculates six types of pair-wise dynamic covariates for a 
#'              network over specified, non-overlapping time intervals. Covariates
#'              are calculated as log(count + 1).
#'
#' @param events A data.frame with three columns: sender, receiver, and time. 
#'               Sender and receiver columns must be integers (e.g., from 1 to n).
#'               The time column should be numeric.
#' @param kz A numeric vector of timestamps defining the end points of the intervals.
#'           For kz = c(t1, t2, t3), intervals will be (0, t1], (t1, t2], (t2, t3].
#' @param n An integer representing the total number of nodes in the network.
#'
#' @return A list with two elements:
#'         - `times`: The original input vector `kz`.
#'         - `covariates`: A list of data.frames. Each data.frame corresponds to a 
#'           time interval and contains only the non-zero pair-wise covariates.
#'           The data.frame has columns: sender, receiver, v1, v2, v3, v4, v5, v6.
#'
#' @details
#' This version builds upon v2 with a key change:
#' 1.  **Count-based Covariates**: Covariates now represent the strength of a
#'     relationship, calculated as `log(count + 1)`, where `count` is the number
#'     of times a specific network pattern occurs (e.g., number of messages from
#'     i to j, or number of 2-step paths from i to j).
#' 2.  It maintains sparse output, interval-based calculation, and numeric node assumptions.
#'
calculate_dynamic_covariates_v3 <- function(events, kz, n) {
  
  # --- 1. Initial Setup ---
  colnames(events) <- c("sender", "receiver", "time")
  time_boundaries <- c(0, kz)
  covariates_list <- vector("list", length(kz))
  
  # --- 2. Loop Through Each Time Interval ---
  for (k in 1:length(kz)) {
    
    t_start <- time_boundaries[k]
    t_end <- time_boundaries[k+1]
    
    interval_events <- events[events$time > t_start & events$time <= t_end, ]
    
    if (nrow(interval_events) == 0) {
      covariates_list[[k]] <- data.frame(sender=integer(), receiver=integer(),
                                         v1=numeric(), v2=numeric(), v3=numeric(),
                                         v4=numeric(), v5=numeric(), v6=numeric())
      next
    }
    
    # --- 3. Build Adjacency Matrix with Counts ---
    
    # Create a sparse adjacency matrix A where A[i, j] = number of events from i to j
    # The `sparseMatrix` function automatically sums up values for duplicate (i, j) pairs.
    adj_matrix_counts <- sparseMatrix(i = interval_events$sender, 
                                      j = interval_events$receiver,
                                      dims = c(n, n),
                                      x = 1) # x=1 for each event
    
    # --- 4. Calculate All 6 Covariate Count Matrices ---
    
    # v1 (send: i -> j): Direct count of messages.
    v1_counts <- adj_matrix_counts
    
    # v2 (receive: i <- j): Transpose of the count matrix.
    v2_counts <- t(adj_matrix_counts)
    
    # For path-based covariates, matrix multiplication now sums up the number of paths.
    # A %*% A gives the number of paths of length 2.
    A_squared_counts <- adj_matrix_counts %*% adj_matrix_counts
    
    # v3 (2-send: i -> h -> j): Count of 2-step paths from i to j.
    v3_counts <- A_squared_counts
    
    # v4 (2-receive: i <- h <- j): Transpose of the 2-send count matrix.
    v4_counts <- t(v3_counts)
    
    # v5 (sibling: h -> i, h -> j): Count of common senders.
    # The entry (i, j) in t(A) %*% A is sum_h(A[h,i] * A[h,j]).
    v5_counts <- t(adj_matrix_counts) %*% adj_matrix_counts
    
    # v6 (cosibling: i -> h, j -> h): Count of common receivers.
    # The entry (i, j) in A %*% t(A) is sum_h(A[i,h] * A[j,h]).
    v6_counts <- adj_matrix_counts %*% t(adj_matrix_counts)
    
    # Set diagonals to 0, as covariates are for pairs (i, j) where i != j
    diag(v1_counts) <- diag(v2_counts) <- diag(v3_counts) <- 0
    diag(v4_counts) <- diag(v5_counts) <- diag(v6_counts) <- 0
    
    # --- 5. Generate Sparse Output with Log-Transformed Values ---
    
    # Create a "master mask" matrix to find all pairs (i, j) that have at least one non-zero count
    master_mask <- (v1_counts + v2_counts + v3_counts + v4_counts + v5_counts + v6_counts) > 0
    
    active_pairs_indices <- which(master_mask, arr.ind = TRUE)
    
    if (nrow(active_pairs_indices) == 0) {
      covariates_list[[k]] <- data.frame(sender=integer(), receiver=integer(),
                                         v1=numeric(), v2=numeric(), v3=numeric(),
                                         v4=numeric(), v5=numeric(), v6=numeric())
      next
    }
    
    # Build the sparse data.frame, applying log1p transformation on the fly.
    # log1p(x) is a numerically stable way to compute log(x + 1).
    z_df <- data.frame(
      sender = active_pairs_indices[, "row"],
      receiver = active_pairs_indices[, "col"],
      v1 = log1p(v1_counts[active_pairs_indices]),
      v2 = log1p(v2_counts[active_pairs_indices]),
      v3 = log1p(v3_counts[active_pairs_indices]),
      v4 = log1p(v4_counts[active_pairs_indices]),
      v5 = log1p(v5_counts[active_pairs_indices]),
      v6 = log1p(v6_counts[active_pairs_indices])
    )
    # z_df <- data.frame(
    #   sender = active_pairs_indices[, "row"],
    #   receiver = active_pairs_indices[, "col"],
    #   v1 = log1p(v3_counts[active_pairs_indices]),
    #   v2 = log1p(v4_counts[active_pairs_indices]),
    #   v3 = log1p(v5_counts[active_pairs_indices]),
    #   v4 = log1p(v6_counts[active_pairs_indices])
    # )
    # 
    covariates_list[[k]] <- z_df
  }
  
  # --- 6. Return Final List ---
  
  final_output <- list(
    times = kz,
    covariates = covariates_list
  )
  
  return(final_output)
}

# 确保 Matrix 包已加载
if (!require("Matrix")) {
  install.packages("Matrix")
}
library(Matrix)

#' @title Calculate Dynamic Pair-wise Network Covariates (v2)
#' @description Efficiently calculates six types of pair-wise dynamic covariates for a
#'              network over specified, non-overlapping time intervals.
#'
#' @param events A data.frame with three columns: sender, receiver, and time.
#'               Sender and receiver columns must be integers (e.g., from 1 to n).
#'               The time column should be numeric.
#' @param kz A numeric vector of timestamps defining the end points of the intervals.
#'           For kz = c(t1, t2, t3), intervals will be (0, t1], (t1, t2], (t2, t3].
#' @param n An integer representing the total number of nodes in the network.
#'
#' @return A list with two elements:
#'         - `times`: The original input vector `kz`.
#'         - `covariates`: A list of data.frames. Each data.frame corresponds to a
#'           time interval and contains only the non-zero pair-wise covariates.
#'           The data.frame has columns: sender, receiver, v1, v2, v3, v4, v5, v6.
#'
#' @details
#' This version addresses three key requirements:
#' 1.  **Sparse Output**: Only (sender, receiver) pairs with at least one non-zero
#'     covariate are included in the output data.frames.
#' 2.  **Interval-based Calculation**: Covariates are calculated based on events *within*
#'     each time interval, not cumulatively.
#' 3.  **Numeric Nodes**: Assumes nodes are integers from 1 to n.
#'
#' The term "in the past" from the original image is interpreted as "within the
#' specified time interval" for each calculation.
#'
calculate_dynamic_covariates_v2 <- function(events, kz, n) {

  # --- 1. Initial Setup ---
  colnames(events) <- c("sender", "receiver", "time")

  # Create time boundaries for the intervals. Intervals are (t_start, t_end]
  time_boundaries <- c(0, kz)

  # Prepare the list to store the results for each time interval
  covariates_list <- vector("list", length(kz))


  # --- 2. Loop Through Each Time Interval ---

  for (k in 1:length(kz)) {

    t_start <- time_boundaries[k]
    t_end <- time_boundaries[k+1]

    # Filter events that occurred within the current time interval
    interval_events <- events[events$time > t_start & events$time <= t_end, ]

    # If no events in this interval, create an empty df and move to the next interval
    if (nrow(interval_events) == 0) {
      covariates_list[[k]] <- data.frame(sender=integer(), receiver=integer(),
                                         v1=numeric(), v2=numeric(), v3=numeric(),
                                         v4=numeric(), v5=numeric(), v6=numeric())
      next
    }

    # --- 3. Build Adjacency Matrix for the Interval ---

    # Create a sparse adjacency matrix A for this interval, where A[i, j] = 1 if i -> j
    adj_matrix <- sparseMatrix(i = interval_events$sender,
                               j = interval_events$receiver,
                               dims = c(n, n),
                               x = 1)

    # Ensure the matrix is binary (0 or 1), as multiple events don't change the covariate value
    adj_matrix <- (adj_matrix > 0) * 1

    # --- 4. Calculate All 6 Covariate Matrices ---

    v1_mat <- adj_matrix
    v2_mat <- t(adj_matrix)

    A_squared <- adj_matrix %*% adj_matrix
    v3_mat <- (A_squared > 0) * 1
    v4_mat <- t(v3_mat)

    v5_mat <- (t(adj_matrix) %*% adj_matrix > 0) * 1
    v6_mat <- (adj_matrix %*% t(adj_matrix) > 0) * 1

    # Set diagonals to 0, as covariates are for pairs (i, j) where i != j
    diag(v1_mat) <- diag(v2_mat) <- diag(v3_mat) <- 0
    diag(v4_mat) <- diag(v5_mat) <- diag(v6_mat) <- 0

    # --- 5. Generate Sparse Output ---

    # Create a "master mask" matrix to find all pairs (i, j) that have at least one non-zero covariate
    master_mask <- (v1_mat + v2_mat + v3_mat + v4_mat + v5_mat + v6_mat) > 0

    # Get the (row, col) indices of all non-zero pairs
    active_pairs_indices <- which(master_mask, arr.ind = TRUE)

    # If no pairs have any covariates, create an empty df
    if (nrow(active_pairs_indices) == 0) {
      covariates_list[[k]] <- data.frame(sender=integer(), receiver=integer(),
                                         v1=numeric(), v2=numeric(), v3=numeric(),
                                         v4=numeric(), v5=numeric(), v6=numeric())
      next
    }

    # Efficiently build the sparse data.frame
    z_df <- data.frame(
      sender = active_pairs_indices[, "row"],
      receiver = active_pairs_indices[, "col"],
      v1 = v1_mat[active_pairs_indices],
      v2 = v2_mat[active_pairs_indices],
      v3 = v3_mat[active_pairs_indices],
      v4 = v4_mat[active_pairs_indices],
      v5 = v5_mat[active_pairs_indices],
      v6 = v6_mat[active_pairs_indices]
    )

    covariates_list[[k]] <- z_df
  }

  # --- 6. Return Final List ---

  final_output <- list(
    times = kz,
    covariates = covariates_list
  )

  return(final_output)
}
