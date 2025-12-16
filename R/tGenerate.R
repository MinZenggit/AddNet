#' Generate Simulation Dataset
#'
#' @description Generate the simulation used data using user specific time
#' dependent outgoing and incoming functions. The homofily effect and the
#' non-homofily effect are user provided (currently it only allows non-time
#' dependent effects in this implementation)
#'
#' @param baseline A function take time t as argument and return the intensity
#' at time t.
#' @param fs A function takes two arguments (time, individual id) and return
#' outgoing intensity.
#' @param fr A function takes two arguments (time, individual id) and return
#' incoming intensity.
#' @param fg A function takes two arguments (time, p^th dimension).
#' @param n Number of individuals in the community
#' @param p Number of dimension of non-homofily effect
#' @param zij Non-homofily effect (Need not to be orthogonal to homofily effect)
#' (Only non time-dependent case is supported by this package)
#' @param maxit Maximum individual's intensity (default is 10), user should check
#' carefully because it determines the performance of dataset generation
#'
#' @return
#' \item{trial}{3 column dataframe, first column is the outgoing individual,
#' the second column is the incoming individual while the last column is the
#' time-stamp}
#' \item{n}{Number of individuals in the community}
#' \item{p}{Number of dimension of non-homofily effect}
#' \item{zij}{Non-homofily effect (not the original one) projected using
#' zProjection function internally}
#'
#' @examples
#'
#' # Define Individual's time-dependent intensity
#' n = 50
#' p = 2
#' bs <- function(t) {0.5}
#' fs <- function(t, i) {0.05*(i-(n+1)/2) + 0.05*(i-(n+1)/2)*sin(2*pi*t)}
#' fr <- function(t, i) {0.1*(i-(n+1)/2) + 0.1*(i-(n+1)/2)*sin(2*pi*t)}
#' fg <- function(t, k) {0.2*k}
#'
#' # Constructing Non-homofily effect
#' zij <- array(0, c(n, n, p))
#' for (i in 0:(n-1)) {
#'   for (j in 0:(n-1)) {
#'     if (i < n/3 && j < n/3) {
#'       zij[i, j, 1] = 1
#'     }
#'     if (i >= n/2 && j >= n/2) {
#'       zij[i, j, 2] = 1
#'     }
#'   }
#' }
#'
#' result <- tGenerate(bs, fs, fr, fg, n, p, zij, maxit = 10) # maxit is different by cases
#'
#' @import stats
#' @export
#'
tGenerate <- function(bs, fs, fr, fg, n, p, zij = NULL, tz, maxit = 10) {

  co = 1
  all_sim = list()
  t_l = length(tz)
  t0 = 0
  t1 = tz[1]

  for (iter in seq_len(t_l)) {

    if (!is.null(zij)) {
      if (length(dim(zij[,,,iter])) == 2) {
        zij1 = array(zij[,,,iter], dim = c(dim(zij[,,,iter]), 1))
      } else {
        zij1 = zij[,,,iter]
      }
      zij_t = zProjection(zij1)
    }

    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          kp = rpois(1,maxit*(t1 - t0))
          kt = sort(runif(kp, t0, t1))
          kt_1 = c()
          if (kp != 0) {
            for (z in 1:kp) {
              temp = bs(kt[z]) + fs(kt[z], i) + fr(kt[z], j)
              for (k in seq_len(p)) {
                temp = temp + fg(kt[z], k) * zij_t[i, j, k]
              }
              if (runif(1,0,1) < temp/maxit) {
                kt_1 = c(kt_1, kt[z])
              }
            }
            kt_l = length(kt_1)
            for (z in seq_len(kt_l)) {
              all_sim[[co]] = c(i, j, kt_1[z])
              co = co + 1
            }
          }
        }
      }
    }

    if (iter < t_l) {
      t0 = t1
      t1 = tz[iter + 1]
    }
  }

  trail_sim = data.frame(p = rep(0, length(all_sim)), q = rep(0, length(all_sim)), t2 = rep(0, length(all_sim)))
  for (i in 1:length(all_sim)) {
    trail_sim[i, 1] = all_sim[[i]][1]
    trail_sim[i, 2] = all_sim[[i]][2]
    trail_sim[i, 3] = all_sim[[i]][3]
  }

  rv = order(trail_sim[,3])
  trail_sim = trail_sim[rv,]

  output = list(trail = trail_sim, n = n, p = p, zij = zij)
  output

}
