#' @title Convert data to V statistic
#'
#' @description Convert continuous random variables in in-control process to discrete data with V statistic, where V statistic is the total number of sample satisfying \eqn{Y_{ij}=\frac{(X_{i2j}-X_{i(2j-1)})^2}{2}>\sigma^2} at time \eqn{i},
#' where \eqn{X_{ij}} is the observation for the \eqn{i^{th}} sampling period and the \eqn{j^{th}} sample in the in-control data, \eqn{n} is the number of the sample size and \eqn{m} is the number of the sampling periods.
#' \eqn{\sigma^2} is population variance of continuous in-control data. If \eqn{\sigma^2} is unknown, it can be estimated by \eqn{\hat{\sigma}^2 = \frac{\sum^m_{i=1}S_i^2}{m}} and \eqn{S_i^2 = \frac{\sum^n_{j=1}(X_{ij}-\overline{X}_i)^2}{n-1}}.
#'
#' @param ICdata The in-control data.
#' @param OCdata The out-of-control data.
#' @param var.p Variance of the random variables in the in-control data.
#'
#' @return
#' \code{V0}\eqn{\hspace{2cm}} The V statistic for in-control data.
#'
#' \code{V1}\eqn{\hspace{2cm}} The V statistic for out-of-control data.
#'
#' \code{p0}\eqn{\hspace{2cm}} The process proportion for in-control data.
#'
#' \code{p1}\eqn{\hspace{2cm}} The process proportion for out-of-control data.
#'
#' \code{n}\eqn{\hspace{2.2cm}} The number of the sample size.
#' @export
#'
#' @importFrom stats var
#'
#' @examples
#' IC = matrix(rnorm(100,0,1),ncol = 10,byrow = TRUE)
#' OC = matrix(rnorm(100,0,2),ncol = 10,byrow = TRUE)
#' cont_to_disc_V(IC,OC)
#' @references Yang, S. F. & Arnold, B. C. (2014). A simple approach for monitoring business service time variation.\emph{The Scientific World Journal}, \emph{2014}:16.
#'
#' Yang, S. F., & Arnold, B. C. (2016). A new approach for monitoring process variance. \emph{Journal of Statistical Computation and Simulation}, \emph{86(14)}, 2749-2765.
#'
cont_to_disc_V = function(ICdata,OCdata,var.p = NULL){
  n = ncol(ICdata)
  if(n%%2==1){
    ICdata = ICdata[,-n]
    OCdata = OCdata[,-n]
  }
  if(is.null(var.p)){
    s = mean(apply(ICdata,1,var))
    var.p = s
  }
  n = ncol(ICdata)/2
  m = matrix(ICdata,ncol = 2,byrow = T)
  m1 = matrix(OCdata,ncol = 2,byrow = T)
  Y = (m[,1]-m[,2])^2/2
  Y1 = (m1[,1]-m1[,2])^2/2
  V = matrix(Y>var.p,ncol = n,byrow = T)
  V1 = matrix(Y1>var.p,ncol = n,byrow = T)
  D = apply(V, 1, sum)
  D1 = apply(V1, 1, sum)
  p = mean(apply(V, 1, mean))
  p1 = mean(apply(V1, 1, mean))
  return(list(V0 = D,V1 = D1,
              p0 = p,p1 = p1,n = n))
}
