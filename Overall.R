


#' The p-value of overall method.
#'
#' @param pv the p-values that will be aggregated. In this study, the p-values are the output of 
#'    either Summary_Test or Uncorrect_Test, which represent the p-values of 
#'    all BT, SKAT, SKATO with and without eQTL-derived weights using GWAS summary statistics 
#'    and correct LD structure, the number of output is L = 3*(n_study + 1).
#' @param Omega Estimated correlation matrix of p-values under H0.
#'
#' @return p-value of overall method.



Overall <- function(pv, Omega){
  pv <- as.numeric(pv)
  pv_j <- sort(pv, index.return = TRUE)
  idx <- pv_j$ix
  r <- Omega[idx, idx]
  rho <- r
  m <- length(pv)
  me=numeric(m)
  for(i in 1:m){
    rhoi=rho[1:i,1:i]
    eigeni=eigen(rhoi)$values
    me[i]=i-sum(ifelse(eigeni>1,eigeni,1)-1)
  }   
  pv_overall <- min(me[m]*pv_j$x/me)
  return(pv_overall)
}
