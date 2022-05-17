#' @title A one-sided lower EWMA-p control chart
#'
#' @description This function displays one-sided lower EWMA-p chart control charts based on in-control and out-of-control data that are number of defectives.
#' In the presence of measurement error, this function is able to provide suitable charts with corrections of measurement error effects.
#'
#' @param ICdata The in-control data for attributes.
#' @param OCdata The out-of-control data for attributes.
#' @param lambda An EWMA smooth constant, which is a scalar in [0,1].
#' @param n A sample size in the data.
#' @param pi1 The proportion that the observed defectives are the same as unobserved ones.
#' @param pi2 The proportion that the observed non-defectives are the same as unobserved ones.
#' @param ARL0 A prespecified average run length (ARL) of a control chart in the in-control process.
#' @param M The number of simulation times for the Monte Carlo method
#' @param error The tolerant for the absolute difference between an iterated ARL value and prespecified \code{ARL0}.
#'
#' @return The first chart is an EWMA-p chart obtained by the in-control data, and the second chart is an EWMA-p chart based in the out-of-control data.
#' In two figures, horizontal solid line represents lower control limit (LCL), black solid dots are detections of in-control data, and red solid dots are detections of out-of-control data.
#' @export
#' @importFrom stats rbinom
#' @importFrom stats uniroot
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom graphics abline
#' @importFrom graphics text
#'
#' @references Chen, L.-P. & Yang, S.-F. (2022). A new p-chart with measurement error correction. arXiv: 2203.03384.
#' @examples
#' library(qcr)
#' data = orangejuice
#' IC = data[1:30,1]
#' OC = data[31:54,1]
#' EWMA_p_chart_one_LCL(IC,OC,0.05,50,1,1)
EWMA_p_chart_one_LCL = function(ICdata,OCdata,lambda,n,pi1 = 1,pi2 = pi1,ARL0 = 200,M = 500,error = 10){
  ICdata1 = ICdata/n
  p = mean(ICdata1)
  a = EWMA_p_one_LCL(p,lambda,n,pi1,pi2,ARL0,M,error)
  ICdata1 = (ICdata1+pi2-1)/(pi1+pi2-1)
  E = ewma(ICdata1,lambda,mean(ICdata1))
  color = rep('black',length(ICdata1))
  color[(E-a$LCL)<0] = 'red'
  plot(E,type = 'b',col = color,pch = 16,xlab = 't',ylab = 'EWMA',
       ylim = c(min(c(E,a$LCL))-0.07,max(E)+0.07),main = 'EWMA-p chart for IC data')
  abline(h = a$LCL,lty = 1)
  txt = as.character(round(a$LCL,3))
  text(x = (length(ICdata1)-3), y = a$LCL+0.007, paste('LCL=', txt))

  OCdata1 = OCdata/n
  OCdata1 = (OCdata1+pi2-1)/(pi1+pi2-1)
  E = ewma(OCdata1,lambda,mean(ICdata1))
  color = rep('black',length(OCdata1))
  color[(E-a$LCL)<0] = 'red'
  plot(E,type = 'b',col = color,pch = 16,xlab = 't',ylab = 'EWMA',
       ylim = c(min(c(E,a$LCL))-0.07,max(E)+0.07),main = 'EWMA-p chart for IC data')
  abline(h = a$LCL,lty = 1)
  txt = as.character(round(a$LCL,3))
  text(x = (length(OCdata1)-3), y = a$LCL+0.007, paste('LCL=', txt))
}
