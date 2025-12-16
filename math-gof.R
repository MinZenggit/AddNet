# R anaconda


rm(list = ls())
library(dplyr)
library(ggplot2)
# load(file = "result/re_math.rdata")
load("data/rd_z.rdata")
load("data/rd_y.rdata")
source("codes/format_covariate_tensor.R")
# --- [第1部分：数据准备 - 无变化] ---
# 假设 npara1 和 bike_data 已加载

n <- 1000

# nodes <- sample(1:24818, n, replace = F)
nodes <- 1:1000
format_covariate_tensor_optimized(rd_z$covariates[[1]], nodes) -> zij
trail <- rd_y
names(trail) = c("x", "y", "time")
trail <- trail %>% filter(x %in% nodes, y %in% nodes)
nn = nrow(trail)
p = dim(zij)[3]

source(file = "codes/nonParametric_diff_flexible_2.R")
sourceCpp("codes/create_tensors_flexible_2.cpp")
sourceCpp("codes/compute_flexible_2.cpp")
source(file = "codes/convert_zij_to_sparse_df.R")
z_df = convert_zij_to_sparse_df(zij)
t5 <- Sys.time()
npara1 <- nonParametric_diff_flexible(trail,
                                      z_df,
                                      n, p, K = 100,
                                      smooth = T,
                                      bandwidth = 0.02)
t6 <- Sys.time() 

tseq <- npara1$time_points
alpha_hat <- npara1$truncate$pointed$alpha_hat
beta_hat <- npara1$truncate$pointed$beta_hat
gamma_hat <- npara1$truncate$pointed$gamma_hat
lambda0_hat <- npara1$truncate$pointed$lambda0_hat

# --- [第2部分：创建累积效应函数 - 优化] ---

# 辅助函数保持不变
cumulative_integral <- function(y, x) {
  c(0, cumsum(0.5 * (y[-1] + y[-length(y)]) * diff(x)))
}

# 步骤1: 预计算所有累积效应的数据矩阵
cat("Step 1: Pre-calculating cumulative effect matrices...\n")
Lambda0_val <- cumulative_integral(lambda0_hat, tseq)
# 使用sapply/lapply高效创建矩阵
A_hat_cum <- t(sapply(1:n, function(i) cumulative_integral(alpha_hat[i, ], tseq)))
B_hat_cum <- t(sapply(1:n, function(j) cumulative_integral(beta_hat[j, ], tseq)))
Gamma_hat_cum <- t(sapply(1:p, function(k) cumulative_integral(gamma_hat[k, ], tseq)))

# 步骤2: 预计算累积效应的总和
cat("Step 2: Pre-calculating total cumulative effects...\n")
A_total_cum <- colSums(A_hat_cum)
B_total_cum <- colSums(B_hat_cum)

# 步骤3: 为所有效应和它们的总和创建高效的插值函数
cat("Step 3: Creating interpolation functions...\n")
Lambda0_func <- approxfun(tseq, Lambda0_val, rule = 2)
A_funcs <- lapply(1:n, function(i) approxfun(tseq, A_hat_cum[i, ], rule = 2))
B_funcs <- lapply(1:n, function(j) approxfun(tseq, B_hat_cum[j, ], rule = 2))
Gamma_funcs <- lapply(1:p, function(k) approxfun(tseq, Gamma_hat_cum[k, ], rule = 2))

# **关键加速点**: 为总和创建插值函数
A_total_func <- approxfun(tseq, A_total_cum, rule = 2)
B_total_func <- approxfun(tseq, B_total_cum, rule = 2)

# --- [第3部分：预处理协变量 - 无变化] ---
Z_i_sum <- apply(zij, c(1, 3), sum)
Z_j_sum <- apply(zij, c(2, 3), sum)

# --- [第4部分：发送者效应检验 - 加速] ---
cat("Step 4: Calculating goodness-of-fit for senders (vectorized)...\n")
es_Ni_list <- list()
time_range <- range(tseq)

for (i in 1:n) {
  if (i %% 50 == 0) cat("  Processing sender node:", i, "/", n, "\n")
  
  dd <- trail %>% filter(x == i, time >= time_range[1], time <= time_range[2]) %>% arrange(time)
  if (nrow(dd) == 0) next
  
  dd$N <- 1:nrow(dd)
  ind <- if (nrow(dd) > 200) sample(1:nrow(dd), 200, replace = FALSE) %>% sort() else 1:nrow(dd)
  
  # **向量化**: event_times 是一个向量
  event_times <- dd$time[ind]
  
  # 获取当前节点的函数和Z和
  Ai_func <- A_funcs[[i]]
  Bi_func <- B_funcs[[i]]
  Zi_sum_vec <- Z_i_sum[i, ]
  
  # **向量化计算**: 对整个 event_times 向量进行操作
  Lambda0_vec <- Lambda0_func(event_times)
  Ai_vec <- Ai_func(event_times)
  Bi_vec <- Bi_func(event_times)
  
  # 使用高效的总和函数
  sum_B_vec <- B_total_func(event_times) - Bi_vec
  
  # 计算协变量部分 (矩阵乘法)
  # Gamma_mat 维度: length(event_times) x p
  Gamma_mat <- sapply(Gamma_funcs, function(f) f(event_times))
  # Zi_sum_vec 维度: p x 1
  # 结果 covariate_effect_sum_vec 维度: length(event_times) x 1
  covariate_effect_sum_vec <- Gamma_mat %*% Zi_sum_vec
  
  # **最终向量化公式**
  es_Ni <- (n - 1) * Lambda0_vec + 
    (n - 1) * Ai_vec + 
    sum_B_vec + 
    as.vector(covariate_effect_sum_vec) # 转换为向量
  
  es_Ni_list[[i]] <- data.frame(node = i, N = dd$N[ind], es_Ni = es_Ni, time = event_times)
}

# --- [第5部分：接收者效应检验 - 加速] ---
cat("Step 5: Calculating goodness-of-fit for receivers (vectorized)...\n")
es_Nj_list <- list()
for (j in 1:n) {
  if (j %% 50 == 0) cat("  Processing receiver node:", j, "/", n, "\n")
  
  dd <- trail %>% filter(y == j, time >= time_range[1], time <= time_range[2]) %>% arrange(time)
  if (nrow(dd) == 0) next
  
  dd$N <- 1:nrow(dd)
  ind <- if (nrow(dd) > 200) sample(1:nrow(dd), 200, replace = FALSE) %>% sort() else 1:nrow(dd)
  
  event_times <- dd$time[ind]
  
  Aj_func <- A_funcs[[j]]
  Bj_func <- B_funcs[[j]]
  Zj_sum_vec <- Z_j_sum[j, ]
  
  Lambda0_vec <- Lambda0_func(event_times)
  Aj_vec <- Aj_func(event_times)
  Bj_vec <- Bj_func(event_times)
  
  sum_A_vec <- A_total_func(event_times) - Aj_vec
  
  Gamma_mat <- sapply(Gamma_funcs, function(f) f(event_times))
  covariate_effect_sum_vec <- Gamma_mat %*% Zj_sum_vec
  
  es_Nj <- (n - 1) * Lambda0_vec + 
    sum_A_vec + 
    (n - 1) * Bj_vec + 
    as.vector(covariate_effect_sum_vec)
  
  es_Nj_list[[j]] <- data.frame(node = j, N = dd$N[ind], es_Nj = es_Nj, time = event_times)
}

# --- [第6部分：可视化 - 无变化] ---
cat("Step 6: Generating plots...\n")
all_sender_data <- bind_rows(es_Ni_list)
all_receiver_data <- bind_rows(es_Nj_list)

if (nrow(all_sender_data) > 0) {
  p_sender <- ggplot(all_sender_data, aes(x = es_Ni, y = N)) +
    geom_point(alpha = 0.55, color = "dodgerblue", size = 0.5, shape = 16) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    # geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "navy", se = FALSE, size = 0.5) +
    labs(
      x = "Expected Cumulative Events",
      y = "Observed Cumulative Events"
    ) +
    theme_bw() +xlim(c(0, 500)) +ylim(c(0, 500))
  # coord_equal(
  #   xlim = range(all_sender_data$es_Ni, all_sender_data$N, na.rm = TRUE),
  #   ylim = range(all_sender_data$es_Ni, all_sender_data$N, na.rm = TRUE)
  # )
  print(p_sender)
}
ggsave( "mathplots/gof_sender.png", width = 4, height = 4)


# pdf("mathplots/gof_sender.pdf")
# plot(p_sender)
# dev.off()
if (nrow(all_receiver_data) > 0) {
  p_receiver <- ggplot(all_receiver_data, aes(x = es_Nj, y = N)) +
    geom_point(alpha = 0.55, color = "dodgerblue", size = 0.5, shape = 16) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    # geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "darkred", se = FALSE, size = 0.5) +
    labs(
      x = "Expected Cumulative Events",
      y = "Observed Cumulative Events"
    ) +
    theme_bw() +xlim(c(0, 500)) +ylim(c(0, 500))
  # coord_equal(
  #   xlim = range(all_receiver_data$es_Nj, all_receiver_data$N, na.rm = TRUE),
  #   ylim = range(all_receiver_data$es_Nj, all_receiver_data$N, na.rm = TRUE)
  # )
  print(p_receiver)
}
ggsave("mathplots/gof_receiver.png", width = 4, height = 4)

# pdf("mathplots/gof_receiver.pdf")
# plot(p_receiver)
# dev.off()



# --- [第7部分：边级别的拟合优度检验 (Goodness-of-Fit for Edges)] ---
# cat("\nStep 7: Calculating goodness-of-fit for specific edges...\n")
# 
# =========================================================================
#               Edge Goodness-of-Fit (All Possible Edges)
# =========================================================================
# 
# 步骤 7.1: 从所有可能的边中抽样
# set.seed(1) # 为了结果可复现
# 
# # 创建一个包含所有可能边（不包括自环）的数据框
# all_possible_edges <- expand.grid(x = nodes, y = nodes) %>% 
#   filter(x != y)
# 
# # 确定要抽样的边的数量
# num_to_sample <- min(1000, nrow(all_possible_edges))
# sampled_edges <- all_possible_edges %>% sample_n(num_to_sample)
# 
# cat("  Total possible edges (excluding self-loops):", nrow(all_possible_edges), "\n")
# cat("  Sampling", num_to_sample, "edges for analysis.\n")
# 
# 
# # 步骤 7.2: 遍历抽样边，为有事件和无事件的边分别计算GoF
# es_Nij_list <- list()
# time_range <- range(tseq)
# t_max <- max(time_range) # 整个观测区间的终点
# 
# pb <- txtProgressBar(min = 0, max = nrow(sampled_edges), style = 3) # 添加进度条
# 
# for (k in 1:nrow(sampled_edges)) {
#   
#   i <- sampled_edges$x[k] # 发送者 (sender)
#   j <- sampled_edges$y[k] # 接收者 (receiver)
#   
#   # 筛选出这条特定边上的所有事件
#   dd_edge <- trail %>% 
#     filter(x == i, y == j, time >= time_range[1], time <= time_range[2]) %>% 
#     arrange(time)
#   
#   # --- 获取该边对应的函数和协变量 ---
#   Ai_func <- A_funcs[[i]]
#   Bj_func <- B_funcs[[j]]
#   Zij_vec <- zij[i, j, ] # 协变量向量 (长度为 p)
#   
#   # --- 根据有无事件，采用不同策略 ---
#   if (nrow(dd_edge) > 0) {
#     # CASE 1: 这条边有事件 (与原逻辑相同)
#     
#     # 观察到的累积事件数
#     dd_edge$N <- 1:nrow(dd_edge)
#     event_times <- dd_edge$time
#     
#     # 向量化计算在每个事件发生时的期望事件数 (es_Nij)
#     Lambda0_vec <- Lambda0_func(event_times)
#     Ai_vec      <- Ai_func(event_times)
#     Bj_vec      <- Bj_func(event_times)
#     Gamma_mat   <- sapply(Gamma_funcs, function(f) f(event_times))
#     covariate_effect_vec <- Gamma_mat %*% Zij_vec
#     
#     es_Nij <- Lambda0_vec + Ai_vec + Bj_vec + as.vector(covariate_effect_vec)
#     
#     # 存储结果 (多行)
#     es_Nij_list[[k]] <- data.frame(
#       sender = i, receiver = j, N = dd_edge$N, es_Nij = es_Nij, time = event_times
#     )
#     
#   } else {
#     # CASE 2: 这条边没有事件
#     
#     # 观察到的事件总数 N = 0
#     # 我们需要计算模型在整个时间窗口 [t_min, t_max] 内预测的事件总数
#     
#     Lambda0_total <- Lambda0_func(t_max)
#     Ai_total      <- Ai_func(t_max)
#     Bj_total      <- Bj_func(t_max)
#     Gamma_total_vec <- sapply(Gamma_funcs, function(f) f(t_max))
#     covariate_effect_total <- sum(Gamma_total_vec * Zij_vec)
#     
#     es_Nij_total <- Lambda0_total + Ai_total + Bj_total + covariate_effect_total
#     
#     # 存储结果 (单行)
#     es_Nij_list[[k]] <- data.frame(
#       sender = i, receiver = j, N = 0, es_Nij = es_Nij_total, time = t_max
#     )
#   }
#   setTxtProgressBar(pb, k) # 更新进度条
# }
# close(pb)
# 
# 
# # 步骤 7.3: 合并数据并可视化
# all_edge_data <- bind_rows(es_Nij_list)
# 
# if (nrow(all_edge_data) > 0) {
#   
#   # 动态确定绘图范围，使其为正方形
#   plot_max <- max(all_edge_data$N, all_edge_data$es_Nij, na.rm = TRUE)
#   
#   p_edge <- ggplot(all_edge_data, aes(x = es_Nij, y = N)) +
#     geom_point(alpha = 0.2, color = "darkgreen", size = 1.2, shape = 16) +
#     geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 0.8) +
#     labs(
#       title = "Goodness-of-Fit for Edges (Including Zero-Event Edges)",
#       subtitle = paste("Based on a random sample of", num_to_sample, "from all possible edges"),
#       x = "Expected Cumulative Events on Edge (Model Prediction)",
#       y = "Observed Cumulative Events on Edge (Data)",
#       caption = "Points on the x-axis (y=0) represent edges with zero observed events."
#     ) +
#     theme_bw() +
#     # 使用 coord_fixed 确保 y=x 线是完美的45度角，并动态设置范围
#     coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max))
#   
#   print(p_edge)
#   
#   # 保存PDF
#   ggsave("mathplots/gof_edge_all_possible.pdf", plot = p_edge, width = 7, height = 7)
#   
# } else {
#   cat("No edge data was generated for plotting.\n")
# }
# 
# pdf("mathplots/gof_edge.pdf")
# plot(p_edge)
# dev.off()
# # ========
# # =========================================================================
# # Part 8: 新增功能 - 计算边级别的定量拟合优度指标 (ISE)
# # =========================================================================
# cat("\nStep 8: Calculating quantitative goodness-of-fit metric for edges (ISE)...\n")
# 
# # 检查是否已抽样边，以防出错
# if (exists("sampled_edges") && nrow(sampled_edges) > 0) {
#   
#   # 用于存储每条边的积分平方误差
#   edge_ise_values <- numeric(nrow(sampled_edges))
#   
#   # 创建进度条
#   pb_ise <- txtProgressBar(min = 0, max = nrow(sampled_edges), style = 3)
#   
#   # 遍历之前抽样的200条边
#   for (k in 1:nrow(sampled_edges)) {
#     
#     i <- sampled_edges$x[k]
#     j <- sampled_edges$y[k]
#     
#     # 1. 获取该边的事件数据
#     dd_edge <- trail %>% 
#       filter(x == i, y == j, time >= time_range[1], time <= time_range[2]) %>% 
#       arrange(time)
#     
#     event_times <- dd_edge$time
#     
#     # 2. 定义积分的分割点：包含起始点、所有事件点和结束点
#     integration_points <- unique(sort(c(time_range[1], event_times, time_range[2])))
#     
#     # 3. 为这条边 (i,j) 构建其特定的 N_hat(t) 预测函数
#     #    这个函数会利用之前创建的所有插值函数
#     Ai_func <- A_funcs[[i]]
#     Bj_func <- B_funcs[[j]]
#     Zij_vec <- zij[i, j, ]
#     
#     N_hat_func_ij <- function(t) {
#       # 向量化计算协变量的累积效应
#       # sapply 会对 t 中的每个时间点调用 Gamma_funcs, 返回一个矩阵
#       gamma_vals_at_t <- sapply(Gamma_funcs, function(f) f(t))
#       # 如果 t 是单个值, gamma_vals_at_t 是向量, 结果是标量
#       # 如果 t 是向量, gamma_vals_at_t 是矩阵, 结果是向量
#       covariate_part <- if (is.vector(gamma_vals_at_t)) {
#         sum(gamma_vals_at_t * Zij_vec)
#       } else {
#         gamma_vals_at_t %*% Zij_vec
#       }
#       
#       # 将所有累积效应相加
#       result <- Lambda0_func(t) + Ai_func(t) + Bj_func(t) + as.vector(covariate_part)
#       return(result)
#     }
#     
#     # 初始化这条边的积分误差
#     total_ise_for_edge <- 0
#     
#     # 4. 遍历所有子区间进行数值积分
#     for (m in 1:(length(integration_points) - 1)) {
#       t_start <- integration_points[m]
#       t_end   <- integration_points[m+1]
#       
#       if (t_end <= t_start) next
#       
#       t_mid <- (t_start + t_end) / 2
#       
#       # a) 计算 N_hat(t) 在中点的值
#       N_hat_at_mid <- N_hat_func_ij(t_mid)
#       
#       # b) 计算 N(t) 在该区间内的常数值
#       N_observed_in_interval <- sum(event_times <= t_start)
#       
#       # c) 计算该子区间的积分平方误差并累加
#       squared_error <- (N_hat_at_mid - N_observed_in_interval)^2
#       interval_width <- t_end - t_start
#       total_ise_for_edge <- total_ise_for_edge + squared_error * interval_width
#     }
#     
#     # 存储计算出的该边的总ISE
#     edge_ise_values[k] <- total_ise_for_edge
#     
#     # 更新进度条
#     setTxtProgressBar(pb_ise, k)
#   }
#   
#   close(pb_ise)
#   
#   # 5. 计算所有抽样边的平均积分平方误差
#   mean_ise <- mean(edge_ise_values, na.rm = TRUE)
#   
#   # 6. 打印最终结果
#   cat("\n\n--- Quantitative Goodness-of-Fit Metric ---\n")
#   cat(paste("Metric: Integrated Squared Error (ISE) = ∫(N̂(t) - N(t))² dt\n"))
#   cat(paste("Number of sampled edges:", nrow(sampled_edges), "\n"))
#   cat(sprintf("Average ISE across sampled edges: %.4f\n", mean_ise))
#   cat("-------------------------------------------\n")
#   cat("Note: A lower value indicates a better model fit.\n")
#   
# } else {
#   cat("\nSkipping ISE calculation because no edges were sampled.\n")
# }
