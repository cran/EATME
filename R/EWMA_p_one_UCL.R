#' @title The one-sided upper control limit of an EWMA-p chart
#'
#' @description This function is used to calculate the one-sided upper control limit for EWMA-p charts with the correction of measurement error effects.
#' If two truly classified probabilities \code{pi1} and \code{pi2} are given by 1, then the corresponding control limit is free of measurement error.
#'
#' @param p The proportion of defectives in the in-control process.
#' @param lambda An EWMA smooth constant, which is a scalar in [0,1].
#' @param n A sample size in the data.
#' @param pi1 The proportion that the observed defectives are the same as unobserved ones.
#' @param pi2 The proportion that the observed non-defectives are the same as unobserved ones.
#' @param ARL0 A prespecified average run length (ARL) of a control chart in the in-control process.
#' @param M The number of simulation times for the Monte Carlo method
#' @param error The tolerant for the absolute different between an itevated ARL calue and prespecified \code{ARL0}.
#'
#' @return
#' \code{L1}\eqn{\hspace{2.2cm}} The coefficient of the upper control limit.
#'
#' \code{hat_ARL0}\eqn{\hspace{1.1cm}} The estimated in-control average run length based on given \code{L1}.
#'
#' \code{hat_MRL0}\eqn{\hspace{1.1cm}} The estimated in-control median of run length based on given \code{L1}.
#'
#' \code{hat_SDRL0}\eqn{\hspace{0.9cm}} The estimated in-control standard deviation of run length based on given \code{L1}.
#'
#' \code{UCL}\eqn{\hspace{2cm}} The limiting value of the upper control limit with \code{L1}.
#' @export
#' @importFrom stats rbinom
#' @importFrom stats uniroot
#' @importFrom stats sd
#' @importFrom stats median
#' @references Chen, L. P., & Yang, S. F. (2022). A New \eqn{p}-Control Chart with Measurement Error Correction. \emph{arXiv preprint} arXiv:2203.03384.
#'
#' @examples
#' EWMA_p_one_UCL(0.2,0.05,5,1,1)
EWMA_p_one_UCL = function(p,lambda,n,pi1 = 1,pi2 = pi1,ARL0 = 200,M = 500,error = 10) {
  p0 = p
  p = (p+pi2-1)/(pi1+pi2-1)
  v = p0*(1-p0)/(pi1+pi2-1)^2
  RL = function(lambda,p,L,v,n){
    data = rbinom(1000,n,p)/n
    s = rep(NA,1000)
    t = rep(NA,1000)
    s[1] = lambda * data[1] + (1 - lambda) * p#EWMA_1
    t[1] = p+L*sqrt(lambda*(1-(1-lambda)^2)*v/(2-lambda)/n)#UCL_1
    i = 1
    while(s[i]<t[i]){
      i = i+1
      s[i] = lambda * data[i] + (1 - lambda) * s[i-1];#EWMA_i
      t[i] = p+L*sqrt(lambda*(1-(1-lambda)^(2*(i)))*v/(2-lambda)/n)#EWMA_i
      if (i == length(data)){
        data = c(data,rbinom(1000,n,p)/n)
        s = c(s,rep(NA,1000))
        t = c(t,rep(NA,1000))
      }
    }
    return(i)
  }

  ARLF = function(lambda,p,L,v,n){
    a = replicate(M,RL(lambda,p,L,v,n))
    return(list(ARL = mean(a), MRL = median(a), SDRL = sd(a)))
  }

  af = function(x){
    ARLF(lambda,p,x,v,n)$ARL-ARL0
  }

  L = uniroot(af,lower = 0,upper = 3, extendInt = "yes")
  hat_ARL0 = ARLF(lambda,p,L$root,v,n)
  ARL0_error = abs(hat_ARL0$ARL-ARL0)
  while (ARL0_error>error){
    L = uniroot(af,lower = 0,upper = 3, extendInt = "yes")
    hat_ARL0 = ARLF(lambda,p,L$root,v,n)
    ARL0_error = abs(hat_ARL0$ARL-ARL0)
  }
  return(list(L1 = L$root, hat_ARL0 = hat_ARL0$ARL,
              hat_MRL = hat_ARL0$MRL, hat_SDRL = hat_ARL0$SDRL,
              UCL = p+L$root*sqrt(lambda*v/(2-lambda)/n)))
}
