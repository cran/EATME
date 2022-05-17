#' @title Generate the discrete random variable with measurement error
#'
#' @description Generate the discrete random variable with measurement error.
#'
#' @param p A probability of the unobserved defectives.
#' @param n A number of sample size in the data.
#' @param m A number of observation in the data.
#' @param pi1 The proportion that the observed defectives are the same as unobserved ones.
#' @param pi2 The proportion that the observed non-defectives are the same as unobserved ones.
#'
#' @return
#' \code{real_data}\eqn{\hspace{1.2cm}} The generated data without measurement error.
#'
#' \code{obs_data}\eqn{\hspace{1.4cm}} The generated data with measurement error.
#'
#' \code{n}\eqn{\hspace{2.6cm}} A sample size in the generated data.
#' @export
#'
#' @examples
#' ME_data_generate(0.7,50,50,0.95)
ME_data_generate = function(p,n,m,pi1, pi2 = pi1){
  X = rbinom(n*m,1,p)
  Y = X
  Y_0_index = sample(which(X==1),round((1-pi1)*sum(X==1)))
  Y[Y_0_index] = abs(Y[Y_0_index]-1)
  Y_1_index = sample(which(X==0),round((1-pi2)*sum(X==0)))
  Y[Y_1_index] = abs(Y[Y_1_index]-1)
  X = apply(matrix(X,ncol = n,byrow = T),1,sum)
  Y = apply(matrix(Y,ncol = n,byrow = T),1,sum)
  return(list(real_data=X ,obs_data = Y,n = n))
}
