


#' The p-values of overall method
#'
#' @param pv 
#' @param Omega 
#'
#' @return



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
