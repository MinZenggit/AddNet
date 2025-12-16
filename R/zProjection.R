#' Projection of Non-homogeneous effect
#'
#' @description
#' Construct the non-honogeneous effect on residual space (Proporsition 1)
#'
#'@param A The non-homogeneous effect (array with dimension (n,n,p))
#'
#'@return
#'\item{zij}{The transfromed non-homogeneous effect (array with dimension (n,n,p))}
#'
#'@export
#'
zProjection <- function(A) {
  n = dim(A)[1]
  p = dim(A)[3]
  diff = 3
  while (diff > 1e-6) {
    for (i in seq_len(p)) {
      r1 = rowSums(A[,,i]) / (n-1)
      for (t1 in 1:n) {
        A[t1, , i] = A[t1, , i] - r1[t1]
        A[t1, t1, i] = 0
      }
      r2 = colSums(A[,,i]) / (n-1)
      for (t1 in 1:n) {
        A[, t1, i] = A[, t1, i] - r2[t1]
        A[t1, t1, i] = 0
      }
      diff = sum(abs(r1))
    }
  }
  A
}
