rm(list = ls())
library(matrixcalc)
library(Rcpp)
library(RcppEigen)
library(MASS)
library(ggplot2)
library(latex2exp)
library(dplyr)
library(Matrix)
library(data.table)
# source(file = "nonParamteric.R")
# sourceCpp(file = "src/splinecalc.cpp")
# sourceCpp(file = "src/cumcalc.cpp")
source(file = "R/tGenerateC.R") 
source(file = "R/xConstruct.R")
source(file = "codes/analyze_and_plot.R")
source(file = "codes_tv/nonParametric_diff_flexible_tv.R")
sourceCpp("codes_tv/create_tensors_flexible_tv.cpp")
sourceCpp("codes_tv/compute_flexible_para_tv.cpp")
source("codes_tv/compute_covariates.R")
# source(file = "codes/nonParametric_diff_flexible.R")
# sourceCpp("codes/create_tensors_flexible.cpp")
# sourceCpp("codes/compute_flexible.cpp")
# source(file = "codes/convert_zij_to_sparse_df.R")

fread("~/Additive_network/data/sx-mathoverflow.txt") -> math
math[order(math[, 3]), ] -> math
names(math) <- c("s", "r", "t")
math$s %>% unique() %>% length() # 19774
math$r %>% unique() %>% length() # 21980
max(math$s)
max(math$r)


# 1. 找出所有唯一的节点ID
# 将s列和r列合并，然后取唯一值，并排序（排序是为了让映射结果稳定）
all_nodes <- sort(unique(c(math$s, math$r)))

node_map <- setNames(1:(length(all_nodes)), all_nodes)
math_cleaned <- data.frame(
  s = unname(node_map[as.character(math$s)]),
  r = unname(node_map[as.character(math$r)]),
  t = math$t
)

# 4. 验证结果
cat("\n--- 清理后的结果 ---\n")
cat("新数据中的唯一节点数:", length(unique(c(math_cleaned$s, math_cleaned$r))), "\n")
cat("新数据中的最小节点ID:", min(math_cleaned$s, math_cleaned$r), "\n")
cat("新数据中的最大节点ID:", max(math_cleaned$s, math_cleaned$r), "\n")

head(math_cleaned)
# 1254192988 --> 1457262355
math_cleaned %>% filter(s <= 24818, s > 0, r <= 24818, r > 0, t <= 1457262355, t > 1254192988) -> rd
rd$s %>% unique() -> ss
rd$r %>% unique() -> rr
max(c(ss, rr)) -> n
rd$t <- (rd$t - min(rd$t))/(max(rd$t) - min(rd$t)) * 1.1
time_points = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
rd_z <- calculate_dynamic_covariates_v3(rd,
                                        time_points,
                                        n = n)
rd_z$times <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
rd %>% filter(t > 0.1) -> rd_y
rd_y$t - 0.1 -> rd_y$t

save(rd_y, file = "data/rd_y.rdata")
save(rd_z, file = "data/rd_z.rdata")
# rd$t <- (rd$t - min(rd$t))/(max(rd$t) - min(rd$t)) * 1.05
# time_points = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
#                 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
# rd_z <- calculate_dynamic_covariates_v3(rd,
#                                         time_points, 
#                                         n = n)
# rd_z$times <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
#                 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95)
# rd %>% filter(t > 0.05) -> rd_y
# rd_y$t - 0.05 -> rd_y$t


t1 <- Sys.time()
re_math <- nonParametric_diff_flexible_tv(
  data = rd_y,
  z_list = rd_z,
  n = n,
  p = 6,
  K = 100,
  smooth = T,
  h1 = 0.05, h2 = 0.05
)
t2 <- Sys.time()
save(re_math, file = "result/re_math.rdata")

cat(as.numeric(t2-t1))#1.03 hours


es <- re_math$smooth$pointed$alpha_hat[1, ]
es_sd <- sqrt(re_math$smooth$pointed$alpha_var[1,])

es <- re_math$truncate$integrated$gamma_hat[2, ]
es_sd <- sqrt(re_math$truncate$integrated$gamma_var[2,2,])

es <- re_math$smooth$pointed$gamma_hat[6, ]
es_sd <- sqrt(re_math$smooth$pointed$gamma_var[6,6,])
# 
es <- re_math$truncate$integrated$gamma_hat[6, ]
es_sd <- sqrt(re_math$truncate$integrated$gamma_var[6,6,])
floordf = data.frame(x = re_math$time_points,
                     y = es,
                     yl = es - 1.96 * es_sd,
                     yu = es + 1.96 * es_sd)
key_dates = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
ggplot(floordf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  geom_vline(
    xintercept = key_dates, 
    linetype = "dotted", 
    color = "black", 
    alpha = 0.7
  )

es <- re_math$truncate$integrated$gamma_hat[5, ]
es_sd <- sqrt(re_math$truncate$integrated$gamma_var[5,5,])
re_ = data.frame(x = re_math$time_points,
                 y = es,
                 yl = es - 1.96 * es_sd,
                 yu = es + 1.96 * es_sd)
ggplot(re_, aes(x = x)) +
  geom_step(aes(y = y), color = "red") +  # 主曲线，黑色
  geom_step(aes(y = yl), color = "blue", linetype = "dashed") +  # 下置信区间，蓝色虚线
  geom_step(aes(y = yu), color = "blue", linetype = "dashed") +  # 上置信区间，红色虚线
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, fill = "grey") +  # 置信区间填充，灰色半透明
  geom_vline(xintercept = key_dates, linetype = "dotted", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(x = "Time", y = "xx")

