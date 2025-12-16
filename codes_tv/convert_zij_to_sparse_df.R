# 首先，确保你安装了 tidyr 包
# install.packages("tidyr")
library(tidyr)

#' 将稠密的协变量数组转换为稀疏数据帧
#'
#' @param zij 一个 n x n x p 的三维数组，代表协变量。
#'
#' @return 一个稀疏格式的 data.frame，包含 sender, receiver, 和 p 个协变量列。
#'         只包含至少有一个协变量不为零的节点对。
convert_zij_to_sparse_df <- function(zij) {
  
  # 检查输入的维度
  dims <- dim(zij)
  if (length(dims) != 3) {
    stop("输入 'zij' 必须是一个 n x n x p 的三维数组。")
  }
  n <- dims[1]
  p <- dims[3]
  
  # --- 核心步骤 ---
  # 1. 使用 which() 找到所有非零元素的位置（行、列、协变量索引）
  #    这是最高效的方法，避免了R中的显式循环。
  #    arr.ind = TRUE 返回一个矩阵，每行是 (row, col, p_idx)
  non_zero_coords <- which(zij != 0, arr.ind = TRUE)
  
  # 如果没有任何非零协变量，返回一个空的数据帧
  if (nrow(non_zero_coords) == 0) {
    cat("警告: 协变量数组 'zij' 中所有值都为零。\n")
    # 创建一个结构正确但没有行的空数据帧
    col_names <- c("sender", "receiver", paste0("p", 1:p))
    return(setNames(data.frame(matrix(ncol = p + 2, nrow = 0)), col_names))
  }
  
  # 2. 获取这些非零元素的值
  non_zero_values <- zij[non_zero_coords]
  
  # 3. 创建一个 "长格式" 的临时数据帧
  long_df <- data.frame(
    sender = non_zero_coords[, 1],
    receiver = non_zero_coords[, 2],
    cov_idx = non_zero_coords[, 3], # 协变量的索引 (1, 2, ..., p)
    value = non_zero_values
  )
  
  # 4. 使用 tidyr::pivot_wider() 将长格式转换为我们需要的 "宽格式"
  #    这是将数据从长变宽的标准、高效操作。
  z_df <- pivot_wider(
    long_df,
    names_from = cov_idx,
    values_from = value,
    names_prefix = "p", # 将列命名为 p1, p2, ...
    values_fill = 0     # 关键！如果一个(i,j)对在p1上有值但在p2上没有，
    # pivot_wider 默认会填NA，我们必须用0填充。
  )
  
  return(z_df)
}
