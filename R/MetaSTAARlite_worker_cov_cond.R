#' Generating summary statistics file for conditional analysis using MetaSTAARlite (the "worker" step)
#'
#' The \code{MetaSTAARlite_worker_cov_cond} function takes in genotype, the genotype
#' of variants to be adjusted for in conditional analysis, the object
#' from fitting the null model, variant information and adjusted variant information (unique identifier)
#' to generate the summary statistics file for the given variant-set, adjusting for a given list of variants.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants. If the input genotype matrix
#' is sparse (e.g. \code{dgCMatrix} format), it is assumed that it has been flipped to represent
#' minor allele coding.
#' @param genotype_adj an n*p_adj genotype matrix (dosage matrix) of the target sequence, where n is
#' the sample size and p_adj is the number of genetic variants to be adjusted for in conditional
#' analysis (or a vector of a single variant with length n if p_adj is 1).
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function
#' in the \code{STAARpipeline} package, or the output from \code{fitNullModel} function in the \code{GENESIS} package
#' and transformed using the \code{\link{genesis2staar_nullmodel}} function in the \code{STAARpipeline} package.
#' @param variant_info a data frame of variant information (unique identifier)
#' with p rows (listed in the same order as the columns of \code{genotype}) and should contain
#' the following 4 columns: chromosome (chr), position (pos), reference allele (ref), and alternative allele (alt).
#' @param variant_adj_info a data frame of adjusted variant information (unique identifier)
#' with p_adj rows (listed in the same order as the rows of \code{genotype_adj}) and should contain
#' the following 4 columns: chromosome (chr), position (pos), reference allele (ref), and alternative allele (alt).
#' @return a list with the following members:
#' @return \code{GTPG_cond}: the covariance matrix between all variants in the variant-set (rows)
#' and all variants in the conditional variant-set (columns) (the covariance file
#' for conditional analysis).
#' @return \code{variant_info}: the data frame of variant information (unique identifier)
#' with p rows (listed in the same order as the rows of \code{GTPG_cond}) and 4 columns: chromosome (chr),
#' position (pos), reference allele (ref), and alternative allele (alt).
#' @return \code{variant_adj_info}: the data frame of adjusted variant information (unique identifier)
#' with p_adj rows (listed in the same order as the columns of \code{GTPG_cond}) and 4 columns: chromosome (chr),
#' position (pos), reference allele (ref), alternative allele (alt), score statistic (U), and variance (V).
#' @references Li, X., et al. (2023). Powerful, scalable and resource-efficient
#' meta-analysis of rare variant associations in large whole genome sequencing studies.
#' \emph{Nature Genetics}, \emph{55}(1), 154-164.
#' (\href{https://doi.org/10.1038/s41588-022-01225-6}{pub})
#' @export

MetaSTAARlite_worker_cov_cond <- function(genotype,genotype_adj,obj_nullmodel,variant_info,variant_adj_info){

  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }
  
  variant_info <- as.data.frame(variant_info)
  variant_adj_info <- as.data.frame(variant_adj_info)

  if(dim(genotype)[2] != dim(variant_info)[1]){
    stop(paste0("Dimensions don't match for genotype and variant_info!"))
  }

  if(inherits(genotype_adj, "numeric") || inherits(genotype_adj, "integer")){
    genotype_adj <- matrix(genotype_adj, ncol=1)
  }

  if(inherits(genotype_adj, "sparseMatrix")){
    # geno_missing_imputation: "minor"
    if(any(is.na(genotype_adj@x)))
    {
      AF_adj <- Matrix::colMeans(genotype_adj,na.rm = TRUE)/2
      impute_vectors_adj <- ifelse(AF_adj<=0.5,0,2)
      genotype_adj <- na.replace.sp(genotype_adj,m=impute_vectors_adj)
    }
    genotype_adj <- as.matrix(genotype_adj)
  } else {
    genotype_adj <- matrix_impute(genotype_adj)
  }

  if(dim(genotype_adj)[2] != dim(variant_adj_info)[1]){
    stop(paste0("Dimensions don't match for genotype_adj and variant_adj_info!"))
  }

  if(dim(genotype)[1] != dim(genotype_adj)[1]){
    stop(paste0("Dimensions don't match for genotype and genotype_adj!"))
  }

  if(!inherits(genotype, "sparseMatrix")){
    genotype <- matrix_flip(genotype)$Geno
    genotype <- as(genotype,"dgCMatrix")
  } else{
    # geno_missing_imputation: "minor"
    if(any(is.na(genotype@x)))
    {
      genotype <- na.replace.sp(genotype,is_NA_to_Zero=TRUE)
    }
  }

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

  GTSinvX <- as(t(genotype) %*% obj_nullmodel$Sigma_iX,"matrix")
  G_adjTSinvX <- as(t(genotype_adj) %*% obj_nullmodel$Sigma_iX,"matrix")
  GTPG_cond <- as(t(obj_nullmodel$Sigma_i %*% genotype) %*% genotype_adj,"matrix") - GTSinvX %*% obj_nullmodel$cov %*% t(G_adjTSinvX)
  G_condTPG_cond <- as(t(obj_nullmodel$Sigma_i %*% genotype_adj) %*% genotype_adj,"matrix") - G_adjTSinvX %*% obj_nullmodel$cov %*% t(G_adjTSinvX)
  V <- diag(t(genotype_adj) %*% obj_nullmodel$Sigma_i %*% genotype_adj - G_adjTSinvX %*% obj_nullmodel$cov %*% t(G_adjTSinvX))
  U <- as.vector(t(genotype_adj) %*% obj_nullmodel$scaled.residuals) # U is for alt allele count
  rm(genotype,genotype_adj)
  gc()

  return(list(GTPG_cond=GTPG_cond,
              G_condTPG_cond=G_condTPG_cond,
              variant_info=variant_info[,c("chr","pos","ref","alt")],
              variant_adj_info=cbind(variant_adj_info[,c("chr","pos","ref","alt")],U,V)))
}
