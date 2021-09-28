


#' P-values of BT, SKAT, SKATO with and without eQTL-derived weights using GWAS summary statistics 
#'    and uncorrect LD structure that is directly estimated by 1KGP reference panel.
#'    
#'    
#'#' This function depends on "mkatr" packages. installation within R: devtools::install_github("baolinwu/mkatr")
#'    (Reference: https://github.com/baolinwu/mkatr).
#'
#'
#' @param Z the vector of Z scores from GWAS summary statistics
#' @param W_hat n_snp * n_study dimensional matrix of eQTL-derived weights.
#' @param R estimated LD structure of a gene from 1KGP reference panel
#'
#' @return he p-values of all BT, SKAT, SKATO with and without eQTL-derived weights using GWAS summary statistics 
#'    and uncorrect LD structure, the number of output is L = 3*(n_study + 1).


library(mkatr)
Uncorrect_Test <- function(Z, W_hat, R){
  unweight <- sats(Z,R)$p.value
  n_study <- ncol(W_hat)
  p.weight <- list()
  for(i.weight in 1:n_study){
    p.weight[[i.weight]] <- sats(Z,R,W_hat[,i.weight])$p.value
  }
  res <- c(unweight, as.vector(do.call(cbind, p.weight)))
  return(res)
}
