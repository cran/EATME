#' @title EWMA chart statistics of the data
#'
#' @description A conventional exponential weighted moving average (EWMA) charting statistic evaluated by the data.
#'
#' @param data An one-dimensional random variable.
#' @param lambda An EWMA smooth constant, which is a scalar in [0,1].
#' @param EWMA0 A starting point of EWMA charting statistic.
#'
#' @return A vector of EWMA charting statistics of \code{data} at different t times.
#' @export
#'
#' @examples
#' x = rnorm(20,0,1)
#' ewma(x,0.05,0)
ewma = function(data,lambda,EWMA0){
  EWMA = rep(NA,length(data))
  EWMA[1] = lambda * data[1] + (1 - lambda) * EWMA0
  for (i in 2:length(data)) {
    EWMA[i] = lambda * data[i] + (1 - lambda) * EWMA[i-1]
  }
  return(EWMA)
}
