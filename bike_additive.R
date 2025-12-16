rm(list = ls())
library(matrixcalc)
library(Rcpp)
library(RcppEigen)
library(MASS)
library(ggplot2)
library(latex2exp)
library(dplyr)
library(Matrix)
# source(file = "nonParamteric.R")
# sourceCpp(file = "src/splinecalc.cpp")
# sourceCpp(file = "src/cumcalc.cpp")
source(file = "R/tGenerateC.R") 
source(file = "R/xConstruct.R")
source(file = "codes/analyze_and_plot.R")
# source(file = "sparse_codes3/nonParametric_diff_flexible_2.R")
# sourceCpp("sparse_codes3/create_tensors_flexible_2.cpp")
# sourceCpp("sparse_codes3/compute_flexible_2.cpp")
source(file = "codes/nonParametric_diff_flexible.R")
sourceCpp("codes/create_tensors_flexible.cpp")
sourceCpp("codes/compute_flexible.cpp")
source(file = "codes/convert_zij_to_sparse_df.R")
load("data/bike_data3.rdata")
trail = bike_data$trail
zij = bike_data$zij
nodes <- 1:542
zij <- bike_data$zij[nodes,nodes,c(1:4)]
zij[is.na(zij)] = 0
trail <- bike_data$trail %>% filter(x %in% nodes, y %in% nodes)
trail <- trail[order(trail$time), ]
n = dim(zij)[1]
nn = nrow(trail)
p = dim(zij)[3]
bandwidth = 0.005
z_df = convert_zij_to_sparse_df(zij)

event_tensors <- create_event_tensors_flexible_cpp(
  senders = trail[, 1],
  receivers = trail[, 2],
  event_times = trail[, 3],
  n = n, K = 100,
  bandwidth = if(is.null(bandwidth)) 0.1 else bandwidth,
  create_kernel = T
)
# apply(event_tensors$truncate,  c(3), colMeans) -> bb
# apply(event_tensors$truncate,  c(3), rowMeans) -> aa

t1 <- Sys.time()
npara1 <- nonParametric_diff_flexible(trail,
                                      array(zij, c(n,n,p,1)),
                                      n, p, K = 100,
                                      smooth = T,
                                      bandwidth = 0.02)
t2 <- Sys.time()
save(npara1, file = "npara1.rdata")
# source(file = "sparse_codes3/nonParametric_diff_flexible_2.R")
# sourceCpp("sparse_codes3/create_tensors_flexible_2.cpp")
# sourceCpp("sparse_codes3/compute_flexible_2.cpp")
# npara2 = nonParametric_diff_flexible(trail, z_df, n, p, K = 100,
#                                      smooth = T, bandwidth = 0.03)
# t3 <- Sys.time()
print(t2-t1)
# print(t3-t2)
es <- npara1$kernel$pointed$alpha_hat[1, ]
es_sd <- sqrt(npara1$kernel$pointed$alpha_var[1,])

es <- npara1$truncate$integrated$gamma_hat[2, ]
es_sd <- sqrt(npara1$truncate$integrated$gamma_var[2,2,])

es <- npara1$kernel$pointed$gamma_hat[1, ]
es_sd <- sqrt(npara1$kernel$pointed$gamma_var[1,1,])
# 
# es <- npara1$truncate$integrated$alpha_hat[1, ]
# es_sd <- sqrt(npara1$truncate$integrated$alpha_var[1, ] %>% as.vector())

floordf = data.frame(x = npara1$time_points,
                     y = es,
                     yl = es - 1.96 * es_sd,
                     yu = es + 1.96 * es_sd)
key_dates = c(0.39, 0.68, 0.86)
ggplot(floordf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  geom_vline(
    xintercept = key_dates, 
    linetype = "dotted", 
    color = "black", 
    alpha = 0.7
  )
