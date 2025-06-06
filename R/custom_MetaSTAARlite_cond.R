#' Performs conditional meta-analysis of custom variant-sets (masks) using MetaSTAARlite
#'
#' This function performs meta-analysis to detect conditional associations between a
#' quantitative/dichotomous phenotype and a user-specified custom mask
#' adjusting for set of known variants by using the MetaSTAARlite pipeline.
#' For each custom mask, the conditional MetaSTAAR-O p-value is a p-value from an omnibus test
#' that aggregated conditional SKAT-MS(1,25), SKAT-MS(1,1), Burden-MS(1,25), Burden-MS(1,1), ACAT-V-MS(1,25),
#' and ACAT-V-MS(1,1) together with conditional p-values of each test weighted by each annotation
#' using Cauchy method.
#' @param chr an integer which specifies the chromosome number.
#' @param mask_name a character which specifies the name of the mask to be meta-analyzed using MetaSTAARlite.
#' @param sample.sizes a numeric vector with the length of \code{study.names}
#' indicating the sample size of each study.
#' @param custom_sumstat_mask_list a list containing study-specific summary statistics corresponding to the custom mask.
#' @param custom_cov_mask_list a list containing study-specific sparse weighted covariance matrices corresponding to the custom mask.
#' @param custom_cov_cond_mask_list a list containing study-specific summary statistics and covariance matrices
#' corresponding to the custom mask for variants to be conditioned on.
#' @param cov_maf_cutoff a numeric vector with the length of \code{study.names}
#' indicating the maximum minor allele frequency cutoffs under which the sparse weighted
#' covariance files between variants are stored.
#' @param rare_maf_cutoff a numeric value specifying the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff an integer specifying the cutoff of minimum number of variants of meta-analyzing
#' a given variant-set (default = 2).
#' @param effect.cond a character value indicating the effects of variants to be adjusted for
#' in conditional analysis are "homogeneous" or "heterogeneous" (default = "homogeneous").
#' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label}
#' specified in \code{\link{custom_MetaSTAARlite_worker}} (default = TRUE).
#' @param variant_type a character value specifying the type of variant included in the analysis. Choices include
#'  "SNV", "Indel", or "variant" (default = "SNV").
#' @param Use_annotation_weights a logical value which determines if annotations will be used as weights or not (default = TRUE).
#' @param Annotation_name a character vector of annotation names used in MetaSTAARlite (default = NULL).
#' @param silent a logical value which determines if the report of error messages will be suppressed (default = FALSE).
#' @return a list of data frames containing the conditional MetaSTAAR p-values (including MetaSTAAR-O) corresponding to the custom mask.
#' @references Li, X., et al. (2023). Powerful, scalable and resource-efficient
#' meta-analysis of rare variant associations in large whole genome sequencing studies.
#' \emph{Nature Genetics}, \emph{55}(1), 154-164.
#' (\href{https://doi.org/10.1038/s41588-022-01225-6}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting noncoding
#' rare-variant associations of large-scale whole-genome sequencing studies.
#' \emph{Nature Methods}.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @export

custom_MetaSTAARlite_cond <- function(chr,mask_name,
                                      sample.sizes,custom_sumstat_mask_list,custom_cov_mask_list,custom_cov_cond_mask_list,
                                      cov_maf_cutoff,rare_maf_cutoff=0.01,rv_num_cutoff=2,effect.cond = c("homogeneous","heterogeneous"),
                                      check_qc_label=TRUE,variant_type=c("SNV","Indel","variant"),
                                      Use_annotation_weights=TRUE,Annotation_name=NULL,silent=FALSE){

  ## evaluate choices
  effect.cond <- match.arg(effect.cond)
  variant_type <- match.arg(variant_type)

  results <- c()

  mask_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,custom_sumstat_mask_list,custom_cov_mask_list,custom_cov_cond_mask_list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(mask_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- mask_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(mask_test_merge_cond,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,3)
      results_temp[2] <- chr
      results_temp[1] <- paste0(as.character(mask_name),"_cond")
      results_temp[3] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results <- rbind(results,results_temp)
    }
  }

  if(!is.null(results))
  {
    colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
    colnames(results)[1:3] <- c("Mask name","Chr","#SNV")
    colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  return(results)
}
