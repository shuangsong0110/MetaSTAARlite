noncoding_MetaSTAARlite <- function(chr,gene_name,
                                    sample.sizes,noncoding_sumstat_gene_list,noncoding_cov_gene_list,
                                    cov_maf_cutoff,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                    check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                    Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)

  ########################################
  #   Downstream

  results_downstream <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["downstream"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["downstream"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "downstream"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_downstream <- rbind(results_downstream,results_temp)
    }
  }

  if(!is.null(results_downstream))
  {
    colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
    colnames(results_downstream)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_downstream)[(dim(results_downstream)[2]-1):dim(results_downstream)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ########################################
  #   Upstream

  results_upstream <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["upstream"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["upstream"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "upstream"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_upstream <- rbind(results_upstream,results_temp)
    }
  }

  if(!is.null(results_upstream))
  {
    colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
    colnames(results_upstream)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_upstream)[(dim(results_upstream)[2]-1):dim(results_upstream)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ########################################################
  #                UTR

  results_UTR <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["UTR"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["UTR"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "UTR"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_UTR <- rbind(results_UTR,results_temp)
    }
  }

  if(!is.null(results_UTR))
  {
    colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
    colnames(results_UTR)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_UTR)[(dim(results_UTR)[2]-1):dim(results_UTR)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  #############################################
  #   Promoter-CAGE

  results_promoter_CAGE <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["promoter_CAGE"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["promoter_CAGE"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "promoter_CAGE"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_promoter_CAGE <- rbind(results_promoter_CAGE,results_temp)
    }
  }

  if(!is.null(results_promoter_CAGE))
  {
    colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
    colnames(results_promoter_CAGE)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_promoter_CAGE)[(dim(results_promoter_CAGE)[2]-1):dim(results_promoter_CAGE)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ##################################################
  #       Promoter-DHS

  results_promoter_DHS <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["promoter_DHS"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["promoter_DHS"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "promoter_DHS"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_promoter_DHS <- rbind(results_promoter_DHS,results_temp)
    }
  }

  if(!is.null(results_promoter_DHS))
  {
    colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
    colnames(results_promoter_DHS)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_promoter_DHS)[(dim(results_promoter_DHS)[2]-1):dim(results_promoter_DHS)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ###########################################
  #        Enhancer-CAGE

  results_enhancer_CAGE <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["enhancer_CAGE"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["enhancer_CAGE"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "enhancer_CAGE"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_enhancer_CAGE <- rbind(results_enhancer_CAGE,results_temp)
    }
  }

  if(!is.null(results_enhancer_CAGE))
  {
    colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
    colnames(results_enhancer_CAGE)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_enhancer_CAGE)[(dim(results_enhancer_CAGE)[2]-1):dim(results_enhancer_CAGE)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ##################################################
  #       Enhancer-DHS

  results_enhancer_DHS <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["enhancer_DHS"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["enhancer_DHS"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "enhancer_DHS"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_enhancer_DHS <- rbind(results_enhancer_DHS,results_temp)
    }
  }

  if(!is.null(results_enhancer_DHS))
  {
    colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
    colnames(results_enhancer_DHS)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_enhancer_DHS)[(dim(results_enhancer_DHS)[2]-1):dim(results_enhancer_DHS)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ############################################
  #           results

  results_noncoding <- list(upstream=results_upstream,downstream=results_downstream,UTR=results_UTR,
                            promoter_CAGE=results_promoter_CAGE,promoter_DHS=results_promoter_DHS,
                            enhancer_CAGE=results_enhancer_CAGE,enhancer_DHS=results_enhancer_DHS)

  return(results_noncoding)
}
