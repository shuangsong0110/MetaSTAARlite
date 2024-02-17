MetaSTAARlite_worker_sumstat <- function(genotype,obj_nullmodel,variant_info,qc_label=NULL,
                                         annotation_phred=NULL,for_individual_analysis=FALSE){
  
  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }
  
  if(inherits(genotype, "sparseMatrix")){
    genotype <- as.matrix(genotype)
  }
  
  if(dim(genotype)[2] != dim(variant_info)[1]){
    stop(paste0("Dimensions don't match for genotype and variant_info!"))
  }
  
  if(!is.null(qc_label) && dim(variant_info)[1] != length(qc_label)){
    stop(paste0("Dimensions don't match for variant_info and qc_label!"))
  }
  
  if(!is.null(annotation_phred) && dim(annotation_phred)[1] != dim(variant_info)[1]){
    stop(paste0("Dimensions don't match for annotation_phred and variant_info!"))
  }
  
  N <- dim(genotype)[1]
  alt_AC <- as.integer(colSums(2 - genotype))
  MAC <- as.integer(pmin(alt_AC, 2 * N - alt_AC))
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  variant_label <- as.vector(MAF>0)
  Geno <- genotype$Geno[,variant_label,drop=FALSE]
  Geno <- as(Geno,"dgCMatrix")
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