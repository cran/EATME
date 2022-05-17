#' @title Convert data to M statistic
#'
#' @description Convert continuous random variables in in-control process into discrete random variables with M statistic, where M statistic is the total number of samples satisfying \eqn{X_{ij}>\mu} at time \eqn{i},
#' where \eqn{X_{ij}} is the observation for the \eqn{i^{th}} sampling period and the \eqn{j^{th}} sample in the in-control data, \eqn{n} is the number of the sample size and \eqn{m} is the number of the sampling periods.
#' \eqn{\mu} is the population mean of continuous in-control data. If \eqn{\mu} is unknown, it can be estimated by \eqn{\hat{\mu}=\overline{\overline{x}}=\frac{\sum^m_{i = 1}\sum^n_{j=1} X_{ij}}{n\times m}}.
#'
#' @param ICdata The in-control data.
#' @param OCdata The out-of-control data.
#' @param mu.p Mean of the random variable in the in-control data.
#'
#' @return
#' \code{M0}\eqn{\hspace{2cm}} The M statistic for in-control data.
#'
#' \code{M1}\eqn{\hspace{2cm}} The M statistic for out-of-control data.
#'
#' \code{p0}\eqn{\hspace{2cm}} The process proportion for in-control data.
#'
#' \code{p1}\eqn{\hspace{2cm}} The process proportion for out-of-control data.
#'
#' \code{n}\eqn{\hspace{2.2cm}} The number of the sample size.
#' @export
#'
#' @examples
#' IC = matrix(rnorm(100,0,1),ncol = 10,byrow = TRUE)
#' OC = matrix(rnorm(100,2,1),ncol = 10,byrow = TRUE)
#' cont_to_disc_M(IC,OC)
#' @references Yang, S. F., Lin, J. S., & Cheng, S. W. (2011). A new nonparametric EWMA sign control chart. \emph{Expert Systems with Applications}, \emph{38(5)}, 6239-6243.
#'
#' Yang, S. F. & Arnold, B. C. (2014). A simple approach for monitoring business service time variation.\emph{The Scientific World Journal}, \emph{2014}:16.
#'
#' Yang,  S. F. (2016). An improved distribution-free EWMA mean chart. \emph{Communications in Statistics-Simulation and Computation}, \emph{45(4)}, 1410-1427.
#'
cont_to_disc_M = function(ICdata,OCdata,mu.p = mean(ICdata)){
  M = matrix(ICdata>mu.p,ncol = ncol(ICdata),byrow = T)
  M1 = matrix(OCdata>mu.p,ncol = ncol(OCdata),byrow = T)
  D = apply(M, 1, sum)
  D1 = apply(M1, 1, sum)
  p = mean(apply(M, 1, mean))
  p1 = mean(apply(M1, 1, mean))
  return(list(M0 = D,M1 = D1,
              p0 = p,p1 = p1,n = ncol(ICdata)))
}
