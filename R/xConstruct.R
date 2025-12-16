#' Projection and X creation
#'
#' @description
#' Create X matrix
#'
#'@param n Sample size
#'@param A The non-homogeneous effect (array with dimension (n,n,p))
#'
#'@return
#'\item{intercept}{Intercept}
#'\item{x}{Homofily effect}
#'\item{zij}{The transfromed non-homogeneous effect (array with dimension (n,n,p))}
#'
#'@export
#'

xConstruct <- function(n, A = NULL) {
  x = array(0, dim = c((n*(n-1)), 2*(n - 1)))

  if (!is.null(A)) {
    p = dim(A)[3]
    t = dim(A)[4]
    zij = array(0, c(n*(n-1), p, t))

    for (iter in seq_len(t)) {
      co = 1
      if (length(dim(A[,,,iter])) == 2) {
        A1 = array(A[,,,iter], dim = c(dim(A[,,,iter]), 1))
      } else {
        A1 = A[,,,iter]
      }

      for (i in 1:n) {
        for (j in 1:n) {
          if (i != j) {
            for (q in seq_len(p)) {
              zij[co, q, iter] = A1[i, j, q]
            }
            co = co + 1
          }
        }
      }
    }
  } else {
    zij = NULL
  }

  co = 1
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        if (i == 1) {
          for (k in 1:(n-1)) {
            x[co, k] = -1
          }
        } else {
          x[co, i-1] = 1
        }
        if (j == 1) {
          for (k in 1:(n-1)) {
            x[co, k+(n-1)] = -1
          }
        } else {
          x[co, j-1+(n-1)] = 1
        }
        co = co + 1
      }
    }
  }

  list(intercept = matrix(1, nrow = n*(n-1), ncol = 1), x = x, zij = zij)
}
