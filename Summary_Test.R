








#' P-values of BT, SKAT, SKATO with and without eQTL-derived weights using GWAS summary statistics 
#'    and correct LD structure.
#'    
#' This function depends on "sumFREGAT" packages. Only availible in Linux with "bcftools" 
#'    (Please see https://github.com/samtools/bcftools to install).
#' 
#' @param snp The information of SNPs in a gene. It must contain "chr", "pos", "ID", 
#'    "EA" (risk allele), "P" (p-value), "beta" (effect size)
#' @param pos.maf The index of position, in which snps with MAF<0.05 will be excluded in this analysis. 
#' @param W_hat n_snp * n_study dimensional matrix of eQTL-derived weights.
#'
#' @return The p-values of all BT, SKAT, SKATO with and without eQTL-derived weights using GWAS summary statistics 
#'    and correct LD structure, the number of output is L = 3*(n_study + 1).
#'


library(sumFREGAT)

Summary_Test <- function(snp, pos.maf, W_hat){
  #------ make weight
  reData <- snp[pos.maf,]
  CHROM <- as.matrix(reData$chr)
  POS <- as.matrix(reData$hg38)
  ID <- as.matrix(reData$SNP_id)
  EA <- as.matrix(reData$a1)
  P <- as.matrix(pchisq(Z^2, 1, lower.tail = FALSE))
  BETA <- beta_hat
  Weights <- rep(1, nrow(reData))
  input.file <- cbind(CHROM, POS, ID, EA, P, BETA, Weights)
  colnames(input.file) <- c("CHROM", "POS", "ID", "EA", "P", "BETA", "Weights")
  write.table(input.file, paste("unweighted_he", n, ".txt", sep = ""), quote = FALSE, row.names = FALSE)
  prep.score.files(input.file = paste("unweighted_he", n, ".txt", sep = ""),
                   output.file.prefix = paste("output_unweighted_he", n, ".scores", sep = ""))
  # system("bgzip output_unweighted.scores.vcf") # zip with bgzip
  system(paste("tabix -p vcf -f output_unweighted_he", n, ".scores.vcf.gz", sep = "")) # to create tabix file
  
  n_study <- ncol(W_hat)
  for (i.weight in 1:n_study){
    Weight <- as.matrix(W_hat[,i.weight])
    input.file <- cbind(CHROM, POS, ID, EA, P, BETA, Weight)
    colnames(input.file) <- c("CHROM", "POS", "ID", "EA", "P", "BETA", "Weights")
    write.table(input.file, paste("Weight",i.weight, "_he", n, ".txt", sep = ""), quote = FALSE, row.names = FALSE)
    prep.score.files(input.file = paste("Weight",i.weight,"_he", n, ".txt", sep = ""),
                     output.file.prefix = paste("output_Weight", i.weight, "_he", n, ".scores", sep = ""))
    # system("bgzip output_CMC.scores.vcf") # zip with bgzip
    system(paste("tabix -p vcf -f output_Weight",i.weight, "_he", n, ".scores.vcf.gz", sep = "")) # to create tabix file
  }
  
  # ------ test association
  # ------------------ unweighted
  gene.file <- "examples/gencode.genePred"
  cor.path <- "examples/cor/"
  score.file <- paste("output_unweighted_he", n, ".scores.vcf.gz", sep = "")
  SKAT.unweighted <- SKAT(score.file, gene.file, genes = gene.name, cor.path = cor.path, 
                          user.weights = TRUE, anno.type = '', beta.par = c(1, 1), 
                          write.file = FALSE, quiet = TRUE)$pvalue
  SKATO.unweighted <- SKATO(score.file, gene.file, genes = gene.name, cor.path = cor.path,
                            user.weights = TRUE, anno.type = '', beta.par = c(1,1),
                            write.file = FALSE, quiet = TRUE)$pvalue
  
  BT.unweighted <- BT(score.file, gene.file, genes = gene.name, cor.path = cor.path,
                      user.weights = TRUE, anno.type = '', beta.par = c(1,1),
                      write.file = FALSE, quiet = TRUE)$pvalue
  
  Test_weight <- list()
  for (i.weight in 1:n_study){
    # ------------------ Weight1
    gene.file <- "examples/gencode.genePred"
    cor.path <- "examples/cor/"
    score.file <- paste("output_Weight", i.weight,"_he", n, ".scores.vcf.gz", sep = "")
    SKAT.weight <- SKAT(score.file, gene.file, genes = gene.name, cor.path = cor.path, 
                        user.weights = TRUE, anno.type = '', beta.par = c(1, 1), 
                        write.file = FALSE, quiet = TRUE)$pvalue
    SKATO.weight <- SKATO(score.file, gene.file, genes = gene.name, cor.path = cor.path,
                          user.weights = TRUE, anno.type = '', beta.par = c(1,1),
                          write.file = FALSE, quiet = TRUE)$pvalue
    
    BT.weight <- BT(score.file, gene.file, genes = gene.name, cor.path = cor.path,
                    user.weights = TRUE, anno.type = '', beta.par = c(1,1),
                    write.file = FALSE, quiet = TRUE)$pvalue
    Test_weight[[i.weight]] <- c(SKATO.weight, SKAT.weight, BT.weight)
  }
  
  result <- c(SKATO.unweighted, SKAT.unweighted, BT.unweighted,
              as.vector(do.call(cbind, Test_weight)))
  return(result)
}
