ncRNA_MetaSTAARlite <- function(chr,gene_name,
                                sample.sizes,ncRNA_sumstat_gene_list,ncRNA_cov_gene_list,
                                cov_maf_cutoff,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){
  
  ## evaluate choices
  variant_type <- match.arg(variant_type)
  
  results <- c()
  
  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,ncRNA_sumstat_gene_list,ncRNA_cov_gene_list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred),silent=silent)
    
    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "ncRNA"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
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
  }
  
  return(results)
}