#' @export
tGenerateC <- function(n, shift1 = 1, shift2 = 0, Zij) {
  
  trail_sim = SimSetC(n, shift1, shift2, Zij)
  
  trail_sim = as.data.frame(trail_sim)
  rv = order(trail_sim[,3])
  trail_sim = trail_sim[rv,]
  nn = nrow(trail_sim)
  colnames(trail_sim) <- c("s", "r", "t")
  
  return(list(trail = trail_sim))
}