#' Generating variant-level summary statistics file using MetaSTAARlite (the "worker" step)
#'
#' The \code{MetaSTAARlite_worker_sumstat} function takes in genotype, the object
#' from fitting the null model, and variant information (unique identifier)
#' to generate the variant-level summary statistics file for the given variant-set.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants. If the input genotype matrix
#' is sparse (e.g. \code{dgCMatrix} format), it is assumed that it has been flipped to represent
#' minor allele coding.
#' @param ALT_AF a numeric vector of alternate allele frequencies for each variant in \code{genotype}. 
#' This is required when \code{genotype} is in sparse format (e.g., \code{dgCMatrix}). 
#' This vector is the \code{ALT_AF} column from the \code{results_info} data frame returned by the \code{Genotype_flip_sp_extraction} function (default = NULL).
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function
#' in the \code{\link{STAARpipeline}} package, or the output from \code{fitNullModel} function in the \code{GENESIS} package
#' and transformed using the \code{\link{genesis2staar_nullmodel}} function in the \code{\link{STAARpipeline}} package.
#' @param variant_info a data frame or matrix of variant information (unique identifier)
#' with p rows (listed in the same order as the columns of \code{genotype}) and should contain
#' the following 4 columns: chromosome (chr), position (pos), reference allele (ref), and alternative allele (alt).
#' @param qc_label a vector of quality control status for each variant in \code{variant_info}, where a PASS variant
#' is labeled as "PASS". If \code{qc_label} is NULL, it is assumed that all variants are PASS variants in the study (default = NULL).
#' @param annotation_phred a data frame or matrix of functional annotation data
#' of dimension p*q (or a vector of a single annotation score with length p).
#' Continuous scores should be given in PHRED score scale, where the PHRED score
#' of j-th variant is defined to be -10*log10(rank(-score_j)/total) across the genome. (Binary)
#' categorical scores should be taking values 0 or 1, where 1 is functional and 0 is
#' non-functional. If not provided, STAAR will perform the
#' SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), ACAT-V(1,1)
#' and ACAT-O tests (default = NULL).
#' @param for_individual_analysis a logical value indicating whether it is used for individual (single-variant) meta-analysis
#' (default = FALSE).
#' @return \code{sumstat}: the data frame of all variants in the variant-set or the list of individual variants
#' (the variant-level summary statistics file), including the following information: chromosome (chr), position (pos), reference allele (ref),
#' alternative allele (alt), quality control status (qc_label, optional), alternative allele count (alt_AC), minor allele count (MAC),
#' minor allele frequency (MAF), study sample size (N), score statistic (U), variance (V), variant annotations provided in
#' \code{annotation_phred}, and the low-rank decomposed component of the covariance file.
#' @references Li, X., et al. (2023). Powerful, scalable and resource-efficient
#' meta-analysis of rare variant associations in large whole genome sequencing studies.
#' \emph{Nature Genetics}, \emph{55}(1), 154-164.
#' (\href{https://doi.org/10.1038/s41588-022-01225-6}{pub})
#' @export

MetaSTAARlite_worker_sumstat <- function(genotype,ALT_AF=NULL,obj_nullmodel,variant_info,qc_label=NULL,
                                         annotation_phred=NULL,for_individual_analysis=FALSE){

  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] != dim(variant_info)[1]){
    stop(paste0("Dimensions don't match for genotype and variant_info!"))
  }

  if(!is.null(qc_label) && dim(variant_info)[1] != length(qc_label)){
    stop(paste0("Dimensions don't match for variant_info and qc_label!"))
  }
  
  if(!is.null(ALT_AF) && dim(variant_info)[1] != length(ALT_AF)){
    stop(paste0("Dimensions don't match for variant_info and ALT_AF!"))
  }

  if(!is.null(annotation_phred) && dim(annotation_phred)[1] != dim(variant_info)[1]){
    stop(paste0("Dimensions don't match for annotation_phred and variant_info!"))
  }

  N <- dim(genotype)[1]
  if(!inherits(genotype, "sparseMatrix"))
  {
    # imputation
    genotype <- matrix_impute(genotype)
    alt_AC <- as.integer(colSums(2 - genotype))
    MAC <- as.integer(pmin(alt_AC, 2 * N - alt_AC))
    genotype <- matrix_flip(genotype)
    MAF <- genotype$MAF
    variant_label <- as.vector(MAF>0)
    Geno <- genotype$Geno[,variant_label,drop=FALSE]
    Geno <- as(Geno,"dgCMatrix")
  } else
  {
    # geno_missing_imputation: "minor"
    if(any(is.na(genotype@x)))
    {
      genotype <- na.replace.sp(genotype,is_NA_to_Zero=TRUE)
    }
    MAC <- Matrix::colSums(genotype)
    alt_AC <- ifelse(ALT_AF<0.5,MAC,2*N-MAC)
    MAF <- MAC/(2*N)
    variant_label <- as.vector(MAF>0)
    Geno <- genotype[,variant_label,drop=FALSE]
  }
  rm(genotype)
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

  GTSinvX_cov <- as(t(Geno) %*% obj_nullmodel$Sigma_iX,"matrix") %*% sqrtm(obj_nullmodel$cov)
  #V <- diag(t(Geno) %*% obj_nullmodel$Sigma_i %*% Geno - GTSinvX_cov %*% t(GTSinvX_cov))
  V <- colSums(Geno * (obj_nullmodel$Sigma_i %*% Geno)) - rowSums(GTSinvX_cov^2) # faster for large number of variants
  U <- as.vector(t(Geno) %*% obj_nullmodel$scaled.residuals)
  rm(Geno)
  gc()

  if (!is.null(qc_label)){
    sumstat <- data.frame(variant_info[,c("chr","pos","ref","alt")],qc_label,alt_AC,MAC,MAF,N)
  }else{
    sumstat <- data.frame(variant_info[,c("chr","pos","ref","alt")],alt_AC,MAC,MAF,N)
  }
  if (!is.null(annotation_phred)){
    sumstat <- cbind(sumstat,annotation_phred)
  }
  sumstat <- sumstat[variant_label,]
  if(!for_individual_analysis){
    sumstat <- cbind(sumstat,U,V,GTSinvX_cov)
  }else{
    sumstat <- cbind(sumstat,U,V)
  }

  return(sumstat)
}
