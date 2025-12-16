library(RcppEigen)
library(ggplot2)
library(Rcpp)
library(matrixcalc)
library(dplyr)
rm(list = ls())
sourceCpp(file = "~/DCCOX/DCCOX/cpp/NT_test.cpp")
# sourceCpp(file = "~/DCCOX/DCCOX/cpp/NT_cv.cpp")
sourceCpp(file = "~/DCCOX/DCCOX/cpp/NT_ini.cpp")
sourceCpp(file = "~/DCCOX/DCCOX/cpp/CI.cpp")
# sourceCpp(file = "cpp/CI_alt.cpp")
sourceCpp(file = "~/DCCOX/DCCOX/cpp/NT_homo.cpp")
source("codes/format_covariate_tensor.R")


load(file = "data/rd_y.rdata")
load(file = "data/rd_z.rdata")
nodes <- 1:1000
format_covariate_tensor_optimized(rd_z$covariates[[1]], nodes) -> zij

trail <- rd_y
names(trail) = c("x", "y", "time")
trail <- trail %>% filter(x %in% nodes, y %in% nodes)
n = dim(zij)[1]
nn = nrow(trail)
p = dim(zij)[3]

h1 = 0.15*n^(-0.1)#0.04375991
h2 = 0.4*n^(-0.2)#0.04787323
tseq = seq(0.005,0.995,0.01)
xkk = matrix(rep(0, 2*n+p-1))

cat(h1, h2, "\n")
# xk1 = NewtonMC_ini(as.matrix(trail), zij, rep(0, 2*n + p - 1), 0.1, h1, h2, n, nn, p)
t1 = Sys.time()
for (t in tseq) {
  xk2 = NewtonMC_test(as.matrix(trail), zij, rep(0.1, 2*n + p - 1), t, h1, h2, n, nn, p)
  xkk = cbind(xkk, matrix(xk2))
  cat(t, "\n")
}
xkk = xkk[, -1]
xkkCI = xkk
count = 1
for (t in tseq) {
  xk = xkk[, count]
  test1 = CI(as.matrix(trail), zij, xk, t, h1, h2, n, nn, p) 
  xkkCI[, count] = test1
  count = count + 1
}
t2 = Sys.time()
print(paste0("n = ", n, ":"))
print((t2-t1)) 
# n = 100; t2-t1 = 34.65706 secs
# n = 300; t2-t1 = 4.345226 mins
# n = 500; t2-t1 = 
# n = 1000; t2-t1 = -----

source(file = "codes/nonParametric_diff_flexible.R")
sourceCpp("codes/create_tensors_flexible.cpp")
sourceCpp("codes/compute_flexible.cpp")
source(file = "codes/convert_zij_to_sparse_df.R")
t3 <- Sys.time()
npara1 <- nonParametric_diff_flexible(trail,
                                      array(zij, c(n,n,p,1)),
                                      n, p, K = 100,
                                      smooth = T,
                                      bandwidth = 0.02)
t4 <- Sys.time() 
# n = 100; t4-t3 = 0.7459683 secs
# n = 300; t4-t3 = 5.225469 secs
# n = 500; t4-t3 = 22.25114 secs
# n = 1000; t4-t3 = 1.45067 mins
print(t4 - t3)

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
print(t6 - t5)
# n = 100; t6-t5 = 0.1879449 secs
# n = 300; t6-t5 = 0.7118213 secs
# n = 500; t6-t5 = 1.867216 secs
# n = 1000; t6-t5 = 7.140288 secs
