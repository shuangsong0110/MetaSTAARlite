#' Generating sparse weighted covariance file using MetaSTAARlite (the "worker" step)
#'
#' The \code{MetaSTAARlite_worker_cov} function takes in genotype and the object
#' from fitting the null model to generate the sparse weighted
#' covariance file for the given variant-set.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants. If the input genotype matrix
#' is sparse (e.g. \code{dgCMatrix} format), it is assumed that it has been flipped to represent
#' minor allele coding.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function
#' in the \code{\link{STAARpipeline}} package, or the output from \code{fitNullModel} function in the \code{GENESIS} package
#' and transformed using the \code{\link{genesis2staar_nullmodel}} function in the \code{\link{STAARpipeline}} package.
#' @param cov_maf_cutoff a numeric value indicating the maximum minor allele frequency cutoff
#' under which the sparse weighted covariance file between variants is stored (default = 0.05).
#' @param qc_label a vector of quality control status for each variant in \code{genotype}, where a PASS variant
#' is labeled as "PASS". If \code{qc_label} is NULL, it is assumed that all variants are PASS variants in the study (default = NULL).
#' @param signif.digits an integer indicating the number of significant digits to be used
#' for storing the sparse weighted covariance file. If \code{signif.digits} is NULL,
#' it is assumed that no rounding will be performed (default = 3).
#' @return \code{GTSinvG_rare}: the sparse matrix of all variants in the variant-set
#' whose minor allele frequency is below \code{cov_maf_cutoff} (the sparse weighted
#' covariance file).
#' @references Li, X., et al. (2023). Powerful, scalable and resource-efficient
#' meta-analysis of rare variant associations in large whole genome sequencing studies.
#' \emph{Nature Genetics}, \emph{55}(1), 154-164.
#' (\href{https://doi.org/10.1038/s41588-022-01225-6}{pub})
#' @export

MetaSTAARlite_worker_cov <- function(genotype,obj_nullmodel,cov_maf_cutoff=0.05,
                                     qc_label=NULL,signif.digits=3){

  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(cov_maf_cutoff < 0 | cov_maf_cutoff > 0.5){
    stop("cov_maf_cutoff should be a number between 0 and 0.5!")
  }

  if(cov_maf_cutoff == 0.5){
    cov_maf_cutoff <- 0.5 + 1e-16
  }

  if(inherits(genotype, "sparseMatrix")){
    MAF <- colMeans(genotype)/2
    if (!is.null(qc_label)){
      RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0)&(qc_label=="PASS"))
    }else{
      RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0))
    }
    Geno_rare <- genotype[,RV_label,drop=FALSE]
    Geno_rare <- as(Geno_rare,"dgCMatrix")
  }else{
    genotype <- matrix_flip(genotype)
    MAF <- genotype$MAF
    if (!is.null(qc_label)){
      RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0)&(qc_label=="PASS"))
    }else{
      RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0))
    }
    Geno_rare <- genotype$Geno[,RV_label,drop=FALSE]
    Geno_rare <- as(Geno_rare,"dgCMatrix")
  }

  rm(genotype,MAF)
  gc()

  if(obj_nullmodel$relatedness){
    if(!obj_nullmodel$sparse_kins){
      stop(paste0("Please use a sparse genetic relatedness matrix when fitting the null model!"))
    }
  }else{
    if(obj_nullmodel$family[1] == "binomial"){
      obj_nullmodel$Sigma_i <- Diagonal(x = obj_nullmodel$weights)
    }else if(obj_nullmodel$family[1] == "gaussian"){
      obj_nullmodel$Sigma_i <- Diagonal(length(obj_nullmodel$y)) / summary(obj_nullmodel)$dispersion
    }
    obj_nullmodel$Sigma_iX <- obj_nullmodel$Sigma_i %*% model.matrix(obj_nullmodel)
    obj_nullmodel$cov <- solve(t(model.matrix(obj_nullmodel)) %*% obj_nullmodel$Sigma_i %*% model.matrix(obj_nullmodel))
    obj_nullmodel$scaled.residuals <- (obj_nullmodel$y - obj_nullmodel$fitted.values) / summary(obj_nullmodel)$dispersion
  }

  GTSinvG_rare <- t(obj_nullmodel$Sigma_i %*% Geno_rare) %*% Geno_rare
  rm(Geno_rare)
  gc()

  GTSinvG_rare <- as(GTSinvG_rare,"TsparseMatrix")
  remove_ind <- (GTSinvG_rare@j < GTSinvG_rare@i) # Be careful about multi-allelic issue
  GTSinvG_rare@i <- GTSinvG_rare@i[!remove_ind]
  GTSinvG_rare@j <- GTSinvG_rare@j[!remove_ind]
  GTSinvG_rare@x <- GTSinvG_rare@x[!remove_ind]
  rm(remove_ind)
  gc()
  GTSinvG_rare <- as(GTSinvG_rare,"CsparseMatrix")

  ### Create a version of GTSinvG_rare with rounded significant digits
  if(!is.null(signif.digits)){
    GTSinvG_rare <- signif(GTSinvG_rare, digits = signif.digits)
  }

  return(GTSinvG_rare)
}
