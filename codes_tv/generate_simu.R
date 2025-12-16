#' @title 生成N阶段时变协变量的模拟数据（带多轮准备阶段）
#'
#' @description
#' 该函数首先执行多个“准备阶段”来生成初始的协变量，然后基于此进行可配置数量的正式模拟阶段。
#' - **准备阶段 (Warm-up)**: 循环执行 `num_warmup_stages` 次。
#'   - 每次循环都基于前一轮的事件计数更新协变量 Z = log(N_ij + 1)。
#'   - 所有准备阶段的数据都将被丢弃，只保留最后一轮生成的Z矩阵。
#' - **正式模拟阶段**:
#'   - 总共执行 `num_formal_stages` 个阶段，时间从 0 到 1 均匀分割。
#'   - 第一个正式阶段使用准备阶段最终生成的 Z 矩阵。
#'   - 后续每个阶段的 Z 都是基于前一个正式阶段的事件计数动态生成的。
#'
#' @param n 节点数。
#' @param p 协变量维度。**注意：当前实现强制 p=1**。
#' @param Zij_initial 一个 n x n x p 的数组，作为**第一轮准备阶段**的协变量。
#' @param max_lambda 强度函数的上界，用于瘦身算法。
#' @param num_warmup_stages 准备阶段的迭代次数。默认为 3。
#' @param num_formal_stages 正式模拟的阶段数量。默认为 5。
#'
#' @return 一个列表，包含正式模拟阶段的结果：
#'         - `events`: 一个包含所有模拟事件的数据框 (sender, receiver, time)，时间范围为 [0, 1]。
#'         - `z_variants`: 一个 n x n x p x num_formal_stages 的数组。
#'         - `z_times`: 协变量变化的 num_formal_stages 个时间点。
#'
generate_simulation_data <- function(n, p, Zij_initial, max_lambda = 10.0, 
                                     num_warmup_stages = 3, num_formal_stages = 5) {
  
  if (p != 1) {
    stop("This simulation function is designed for p=1, as log(count+1) generates a single covariate.")
  }
  if (any(dim(Zij_initial) != c(n, n, p))) {
    stop("Dimensions of Zij_initial must be n x n x p.")
  }
  if (num_formal_stages < 1) {
    stop("num_formal_stages must be at least 1.")
  }
  
  # --- 辅助函数：根据事件数据帧计算下一个Z矩阵 ---
  calculate_next_Z_from_events <- function(events_df, n, p) {
    count_matrix <- matrix(0, nrow = n, ncol = n)
    if (nrow(events_df) > 0) {
      events_df$sender <- factor(events_df$sender, levels = 1:n)
      events_df$receiver <- factor(events_df$receiver, levels = 1:n)
      counts <- table(events_df$sender, events_df$receiver)
      sender_ids <- as.numeric(rownames(counts))
      receiver_ids <- as.numeric(colnames(counts))
      count_matrix[sender_ids, receiver_ids] <- counts
    }
    return(array(log(count_matrix + 1), dim = c(n, n, p)))
  }
  
  # --- 准备阶段 (Warm-up Phase) ---
  # cat(sprintf("Running %d Warm-up Stage(s)...\n", num_warmup_stages))
  Z_warmup_current <- Zij_initial
  warmup_period_duration <- 1 / num_formal_stages # 保持准备阶段和正式阶段的时间段长度一致
  
  if (num_warmup_stages > 0) {
    for (i in 1:num_warmup_stages) {
      events_warmup_raw <- generate_events_for_period_cpp(
        n = n, 
        t_start = 0, 
        t_end = warmup_period_duration,
        Z_period = Z_warmup_current,
        max_lambda = max_lambda,
        p = p
      )
      events_warmup <- as.data.frame(events_warmup_raw)
      Z_warmup_current <- calculate_next_Z_from_events(events_warmup, n, p)
    }
  }
  
  # 准备阶段结束后，Z_warmup_current 中存储的是最后一轮的结果
  Z_formal_start <- Z_warmup_current
  # cat("Warm-up complete. Starting formal simulation with %d stages.\n", num_formal_stages)
  
  # --- 正式模拟开始 (使用循环) ---
  time_breaks <- (0:num_formal_stages) / num_formal_stages
  all_events_list <- list()
  z_variants_list <- list()
  
  Z_current <- Z_formal_start
  
  for (i in 1:num_formal_stages) {
    # cat(sprintf(" -> Simulating Formal Period %d/%d: t in [%.2f, %.2f]...\n", i, num_formal_stages, time_breaks[i], time_breaks[i+1]))
    
    # 1. 存储当前阶段使用的Z矩阵
    z_variants_list[[i]] <- Z_current
    
    # 2. 生成当前阶段的事件
    current_events_raw <- generate_events_for_period_cpp(
      n = n, 
      t_start = time_breaks[i], 
      t_end = time_breaks[i+1],
      Z_period = Z_current,
      max_lambda = max_lambda,
      p = p
    )
    current_events <- as.data.frame(current_events_raw)
    all_events_list[[i]] <- current_events
    
    # 3. 基于刚生成的事件，计算下一个阶段要用的Z矩阵
    Z_current <- calculate_next_Z_from_events(current_events, n, p)
  }
  
  # --- 组合最终结果 ---
  final_events <- do.call(rbind, all_events_list)
  if (nrow(final_events) > 0) {
    final_events <- final_events[order(final_events$time), ]
    final_events$time_scaled <- final_events$time * 1000 
  } else {
    final_events <- data.frame(sender=integer(), receiver=integer(), time=numeric(), time_scaled=numeric())
  }
  
  final_z_variants <- abind::abind(z_variants_list, along = 4)
  final_z_times <- time_breaks[1:num_formal_stages]
  
  return(list(
    events = final_events,
    z_variants = final_z_variants,
    z_times = final_z_times
  ))
}
