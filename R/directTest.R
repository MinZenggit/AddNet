directTest <- function(inc, incvar, out, outvar, tsep) {
  n = nrow(inc)
  s0 = -100
  for (i in tsep) {
    s = sum((inc[, i] - out[, i])^2/(incvar[, i] + outvar[, i]) - 1) / sqrt(n) / sqrt(2)
    s0 = max(s0, s)
  }
  return(list(tstat = s0,
              ind = s0 > qnorm(0.95^(1/length(tsep)))))
}

degreeHTest <- function(inc, incvar, tsep) {
  n = nrow(inc)
  s0 = -100
  for (i in tsep) {
    s = sum((inc[, i]-mean(inc[, i]))^2/(incvar[, i]) - 1) / sqrt(n) / sqrt(2)
    s0 = max(s0, s)
  }
  return(list(tstat = s0,
              ind = s0 > qnorm(0.95^(1/length(tsep)))))
}
