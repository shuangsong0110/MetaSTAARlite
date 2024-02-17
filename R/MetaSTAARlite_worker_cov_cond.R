MetaSTAARlite_worker_cov_cond <- function(genotype,genotype_adj,obj_nullmodel,variant_info,variant_adj_info){
  
  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }
  
  if(dim(genotype)[2] != dim(variant_info)[1]){
    stop(paste0("Dimensions don't match for genotype and variant_info!"))
  }
  
  if(inherits(genotype_adj, "numeric") || inherits(genotype_adj, "integer")){
    genotype_adj <- matrix(genotype_adj, ncol=1)
  }
  
  if(inherits(genotype_adj, "sparseMatrix")){
    genotype_adj <- as.matrix(genotype_adj)
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