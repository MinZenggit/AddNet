#' 将长格式的网络边数据高效转换为三维协变量张量 (Z) [优化版]
#'
#' @description
#' 该函数接收一个以“边列表”形式存储的网络数据框，并将其高效地转换为
#' 一个 n x n x p 的三维数组（张量）。此版本通过完全向量化的三维索引
#' 替代了 for 循环，以获得极致的性能。
#'
#' @param data 一个数据框或矩阵。前两列必须是 sender 和 receiver 的标识符。
#'   其余 p 列是对应的协变量值。
#' @param nodes 一个向量，包含网络中所有节点的唯一标识符。
#'
#' @return 一个 n x n x p 的三维数组 Z。
#'
#' @examples
#' # 创建一个大规模的示例数据集
#' n_nodes <- 500
#' n_edges <- 50000
#' n_covariates <- 6
#' all_nodes <- 1:n_nodes
#' 
#' large_data <- data.frame(
#'   sender   = sample(all_nodes, n_edges, replace = TRUE),
#'   receiver = sample(all_nodes, n_edges, replace = TRUE)
#' )
#' # 添加协变量
#' cov_matrix <- matrix(rnorm(n_edges * n_covariates), ncol = n_covariates)
#' colnames(cov_matrix) <- paste0("v", 1:n_covariates)
#' large_data <- cbind(large_data, cov_matrix)
#' 
#' # 过滤掉自身环路
#' large_data <- large_data[large_data$sender != large_data$receiver, ]
#'
#' # 使用优化后的函数进行基准测试
#' system.time(
#'   Z_tensor <- format_covariate_tensor_optimized(data = large_data, nodes = all_nodes)
#' )
#' # 查看结果
#' # dim(Z_tensor) # 应为 500 500 6

format_covariate_tensor_optimized <- function(data, nodes) {
  
  # --- 1. 输入验证 (与原版相同) ---
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` 必须是一个数据框或矩阵。")
  }
  if (ncol(data) < 3) {
    stop("`data` 必须至少有3列 (sender, receiver, 1个协变量)。")
  }
  if (!is.vector(nodes) || is.list(nodes)) {
    stop("`nodes` 必须是一个向量。")
  }
  
  # --- 2. 初始化 (与原版相同) ---
  unique_nodes <- sort(unique(nodes))
  n <- length(unique_nodes)
  p <- ncol(data) - 2
  
  covariate_names <- colnames(data)[-c(1, 2)]
  if (is.null(covariate_names)) {
    covariate_names <- paste0("v", 1:p)
  }
  
  Z <- array(
    0,
    dim = c(n, n, p),
    dimnames = list(unique_nodes, unique_nodes, covariate_names)
  )
  
  # --- 3. 数据过滤与准备 (与原版相同) ---
  valid_rows_mask <- (data[[1]] %in% unique_nodes) & (data[[2]] %in% unique_nodes)
  
  if (sum(valid_rows_mask) == 0) {
    warning("数据中没有任何边的 sender 和 receiver 同时存在于 `nodes` 向量中。返回一个零张量。")
    return(Z)
  }
  
  filtered_data <- data[valid_rows_mask, , drop = FALSE]
  
  # --- 4. 完全向量化的高效映射与填充 ---
  
  # 将节点ID转换为数组索引
  i_indices <- match(filtered_data[[1]], unique_nodes)
  j_indices <- match(filtered_data[[2]], unique_nodes)
  
  # 获取有效边的数量
  M <- length(i_indices)
  
  # *** 核心优化点：构建三维索引矩阵，取代 for 循环 ***
  
  # 1. 创建一个 (M * p) x 3 的矩阵用于三维索引
  #    第一列是 i (行)，第二列是 j (列)，第三列是 k (切片/协变量)
  full_3d_idx <- cbind(
    rep(i_indices, times = p),
    rep(j_indices, times = p),
    rep(1:p, each = M)
  )
  
  # 2. 将所有协变量值提取并拉直成一个向量
  #    as.matrix() 确保即使 p=1 也能正常工作
  #    as.vector() 会按列展开矩阵，这与我们上面构造 k 索引的方式完美匹配
  covariate_values <- as.vector(as.matrix(filtered_data[, -(1:2)]))
  
  # 3. 使用三维索引进行一次性、完全向量化的赋值
  Z[full_3d_idx] <- covariate_values
  
  # --- 5. 返回结果 ---
  return(Z)
}
