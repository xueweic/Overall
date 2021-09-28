

#' The p-value of copula method.
#'
#' @param pv the p-values that will be aggregated. In this study, the p-values are the output of 
#'    either Summary_Test or Uncorrect_Test, which represent the p-values of 
#'    all BT, SKAT, SKATO with and without eQTL-derived weights using GWAS summary statistics 
#'    and correct LD structure, the number of output is L = 3*(n_study + 1).
#' @param Omega1 Estimated correlation matrix of q-values (inverse-normal p-values) under H0.
#'
#' @return p-value of copula method.


library(mvtnorm)

Copula <- function(pv, Omega1){
  pv_Copula <- 1- pmvnorm(lower = rep(-Inf, length(pv)),
                           upper = rep(qnorm(min(pv), lower.tail = FALSE), length(pv)),
                           mean = rep(0, length(pv)),
                           sigma = Omega1)[1]
  return(pv_Copula)
}
