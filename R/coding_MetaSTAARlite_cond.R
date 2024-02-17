coding_MetaSTAARlite_cond <- function(chr,gene_name,genes,
                                      sample.sizes,coding_sumstat_gene_list,coding_cov_gene_list,coding_cov_cond_gene_list,
                                      cov_maf_cutoff,rare_maf_cutoff=0.01,rv_num_cutoff=2,effect.cond = c("homogeneous","heterogeneous"),
                                      check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                      Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){
  
  ## evaluate choices
  effect.cond <- match.arg(effect.cond)
  variant_type <- match.arg(variant_type)
  
  ### Gene
  kk <- which(genes[,1]==gene_name)
  
  ################################################
  #                  plof_ds
  ################################################
  results_plof_ds <- c()
  
  sumstat.list <- lapply(coding_sumstat_gene_list, function(x) {
    x[["plof_ds"]]
  })
  cov.list <- lapply(coding_cov_gene_list, function(x) {
    x[["plof_ds"]]
  })
  covcond.list <- lapply(coding_cov_cond_gene_list, function(x) {
    x[["plof_ds"]]
  })
  
  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred),silent=silent)
    
    if(inherits(pvalues,"list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "plof_ds_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant
      
      
      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)
      
      results_plof_ds <- rbind(results_plof_ds,results_temp)
    }
  }
  
  if(!is.null(results_plof_ds))
  {
    colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
    colnames(results_plof_ds)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_plof_ds)[(dim(results_plof_ds)[2]-1):dim(results_plof_ds)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }
  
  #####################################################
  #                      plof
  #####################################################
  results_plof <- c()
  
  sumstat.list <- lapply(coding_sumstat_gene_list, function(x) {
    x[["plof"]]
  })
  cov.list <- lapply(coding_cov_gene_list, function(x) {
    x[["plof"]]
  })
  covcond.list <- lapply(coding_cov_cond_gene_list, function(x) {
    x[["plof"]]
  })
  
  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred),silent=silent)
    
    if(inherits(pvalues,"list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "plof_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant
      
      
      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)
      
      results_plof <- rbind(results_plof,results_temp)
    }
  }
  
  if(!is.null(results_plof))
  {
    colnames(results_plof) <- colnames(results_plof, do.NULL = FALSE, prefix = "col")
    colnames(results_plof)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_plof)[(dim(results_plof)[2]-1):dim(results_plof)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }
  
  #############################################
  #             synonymous
  #############################################
  results_synonymous <- c()
  
  sumstat.list <- lapply(coding_sumstat_gene_list, function(x) {
    x[["synonymous"]]
  })
  cov.list <- lapply(coding_cov_gene_list, function(x) {
    x[["synonymous"]]
  })
  covcond.list <- lapply(coding_cov_cond_gene_list, function(x) {
    x[["synonymous"]]
  })
  
  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred),silent=silent)
    
    if(inherits(pvalues,"list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "synonymous_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant
      
      
      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)
      
      results_synonymous <- rbind(results_synonymous,results_temp)
    }
  }
  
  if(!is.null(results_synonymous))
  {
    colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
    colnames(results_synonymous)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_synonymous)[(dim(results_synonymous)[2]-1):dim(results_synonymous)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }
  
  #################################################
  #        missense
  #################################################
  results <- c()
  
  sumstat.list <- lapply(coding_sumstat_gene_list, function(x) {
    x[["missense"]]
  })
  cov.list <- lapply(coding_cov_gene_list, function(x) {
    x[["missense"]]
  })
  covcond.list <- lapply(coding_cov_cond_gene_list, function(x) {
    x[["missense"]]
  })
  
  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred),silent=silent)
    
    if(inherits(pvalues,"list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "missense_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant
      
      
      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)
      
      results <- rbind(results,results_temp)
    }
  }
  
  #################################################
  #         disruptive missense
  #################################################
  sumstat.list <- lapply(coding_sumstat_gene_list, function(x) {
    x[["disruptive_missense"]]
  })
  cov.list <- lapply(coding_cov_gene_list, function(x) {
    x[["disruptive_missense"]]
  })
  covcond.list <- lapply(coding_cov_cond_gene_list, function(x) {
    x[["disruptive_missense"]]
  })
  
  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred),silent=silent)
    
    if(inherits(pvalues,"list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "disruptive_missense_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant
      
      
      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)
      
      results <- rbind(results,results_temp)
    }
  }
  
  if(!is.null(results))
  {
    colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
    colnames(results)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
    
    if(dim(results)[1]==1)
    {
      if(results[3]!="disruptive_missense_cond")
      {
        results <- cbind(results,matrix(1,1,6))
        colnames(results)[(dim(results)[2]-5):dim(results)[2]] <- c("SKAT-MS(1,25)-Disruptive","SKAT-MS(1,1)-Disruptive","Burden-MS(1,25)-Disruptive","Burden-MS(1,1)-Disruptive","ACAT-V-MS(1,25)-Disruptive","ACAT-V-MS(1,1)-Disruptive")
        results_missense <- results
        results_ds <- c()
      }else
      {
        results_missense <- c()
        results_ds <- results
        results <- c()
      }
    }
    
    if(!is.null(results))
    {
      if(dim(results)[1]==2)
      {
        results_m <- c(results[1,],rep(0,6))
        names(results_m)[(length(results_m)-5):length(results_m)] <- c("SKAT-MS(1,25)-Disruptive","SKAT-MS(1,1)-Disruptive","Burden-MS(1,25)-Disruptive","Burden-MS(1,1)-Disruptive","ACAT-V-MS(1,25)-Disruptive","ACAT-V-MS(1,1)-Disruptive")
        results_m[(length(results_m)-5):length(results_m)] <- results[2,c("SKAT-MS(1,25)","SKAT-MS(1,1)","Burden-MS(1,25)","Burden-MS(1,1)","ACAT-V-MS(1,25)","ACAT-V-MS(1,1)")]
        apc_num <- (length(results_m)-18)/6
        p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),1:apc_num+2*(apc_num+1),1:apc_num+3*(apc_num+1),1:apc_num+4*(apc_num+1),1:apc_num+5*(apc_num+1),(6*apc_num+9):(6*apc_num+14))
        results_m["MetaSTAAR-O"] <- CCT(as.numeric(results_m[5:length(results_m)][p_seq]))
        results_m["MetaSTAAR-S(1,25)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num,6*apc_num+9)]))
        results_m["MetaSTAAR-S(1,1)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
        results_m["MetaSTAAR-B(1,25)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
        results_m["MetaSTAAR-B(1,1)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
        results_m["MetaSTAAR-A(1,25)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
        results_m["MetaSTAAR-A(1,1)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))
        
        results_ds <- c()
        results_ds <- rbind(results_ds,results[2,])
        
        results <- c()
        results <- rbind(results,results_m)
      }
    }
  }else
  {
    results <- c()
    results_ds <- c()
  }
  
  results_coding <- list(plof=results_plof,plof_ds=results_plof_ds,missense=results,disruptive_missense=results_ds,synonymous=results_synonymous)
  
  return(results_coding)
}