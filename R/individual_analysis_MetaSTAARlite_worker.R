#' Generates summary statistics of individual variants using MetaSTAARlite
#'
#' This function uses MetaSTAARlite to generate variant-level summary statistics
#' for individual variants of interest.
#' @param chr an integer which specifies the chromosome number.
#' @param start_loc an integer which specifies the starting location of variants to be analyzed.
#' @param end_loc an integer which specifies the end location of variants to be analyzed.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function
#' in the \code{\link{STAARpipeline}} package, or the output from \code{fitNullModel} function in the \code{GENESIS} package
#' and transformed using the \code{\link{genesis2staar_nullmodel}} function in the \code{\link{STAARpipeline}} package.
#' @param known_loci the data frame of variants to be adjusted for in conditional analysis. Should
#' contain four columns in the following order: chromosome (CHR), position (POS), reference allele (REF),
#' and alternative allele (ALT) (default = NULL).
#' @param subsegment.size a numeric value which specifies the size of each subsegment for computation (default = 5e4).
#' @param QC_label a character specifying the channel name of the QC label in the GDS/aGDS file
#' (default = "annotation/filter").
#' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label}.
#' If \code{check_qc_label} is FALSE, then the summary statistics will be stored for PASS variants from the study.
#' If \code{check_qc_label} is FALSE, then the summary statistics will be stored for all variants from the study,
#' together will an additional column of \code{qc_label} (default = FALSE).
#' @param variant_type a character value specifying the type of variant included in the analysis. Choices include
#'  "SNV", "Indel", or "variant" (default = "variant").
#' @param Annotation_dir a character specifying the channel name of the annotations in the aGDS file
#' (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the annotation name and the corresponding
#' channel name in the aGDS file.
#' @param Use_annotation_weights a logical value which specifies if annotations will be used as weights
#' or not (default = FALSE).
#' @param Annotation_name a character vector of annotation names used in MetaSTAARlite (default = NULL).
#' @param silent a logical value which determines if the report of error messages will be suppressed (default = FALSE).
#' @return a list of the following objects corresponding to each individual variant in the given genetic region:
#' (1) the data frame of all variants in the variant-set (the variant-level summary statistics file),
#' including the following information: chromosome (chr), position (pos), reference allele (ref),
#' alternative allele (alt), quality control status (qc_label, optional), alternative allele count (alt_AC), minor allele count (MAC),
#' minor allele frequency (MAF), study sample size (N), score statistic (U), variance (V), variant annotations specified in
#' \code{Annotation_name} (optional), and the low-rank decomposed component of the covariance file;
#' (2) the summary statistics and covariance matrices corresponding to the specified gene for variants to be conditioned on
#' in \code{known_loci}.
#' @references Li, X., et al. (2023). Powerful, scalable and resource-efficient
#' meta-analysis of rare variant associations in large whole genome sequencing studies.
#' \emph{Nature Genetics}, \emph{55}(1), 154-164.
#' (\href{https://doi.org/10.1038/s41588-022-01225-6}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting noncoding
#' rare-variant associations of large-scale whole-genome sequencing studies.
#' \emph{Nature Methods}.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @export

individual_analysis_MetaSTAARlite_worker <- function(chr,start_loc,end_loc,genofile,obj_nullmodel,known_loci=NULL,subsegment.size=5e4,
                                                     QC_label="annotation/filter",check_qc_label=FALSE,variant_type=c("variant","SNV","Indel"),
                                                     Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                                     Use_annotation_weights=c(FALSE,TRUE),Annotation_name=NULL,
                                                     silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)

  phenotype.id <- as.character(obj_nullmodel$id_include)

  filter <- seqGetData(genofile, QC_label)
  if(variant_type=="variant")
  {
    if(check_qc_label)
    {
      SNVlist <- TRUE
    }else
    {
      SNVlist <- filter == "PASS"
    }
  }

  if(variant_type=="SNV")
  {
    if(check_qc_label)
    {
      SNVlist <- isSNV(genofile)
    }else
    {
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
  }

  if(variant_type=="Indel")
  {
    if(check_qc_label)
    {
      SNVlist <- !isSNV(genofile)
    }else
    {
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
  }

  position <- as.integer(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")

  if(!is.null(known_loci))
  {
    allele <- seqGetData(genofile, "allele")

    loc_SNV <- c()
    for (i in 1:dim(known_loci)[1]){
      loc_SNV <- c(loc_SNV,which((position==known_loci$POS[i])&(allele==paste0(known_loci$REF[i],",",known_loci$ALT[i]))))
    }

    seqSetFilter(genofile,variant.id=variant.id[loc_SNV],sample.id=phenotype.id)

    G_SNV <- seqGetData(genofile, "$dosage")

    if (!is.null(G_SNV)){
      id.SNV <- seqGetData(genofile,"sample.id")
      id.SNV.match <- rep(0,length(id.SNV))

      for(i in 1:length(id.SNV))
      {
        id.SNV.match[i] <- which.max(id.SNV==phenotype.id[i])
      }

      G_SNV <- G_SNV[id.SNV.match,,drop=FALSE]
      G_SNV_MAF <- apply(G_SNV, 2, function(x){sum(x[!is.na(x)])/2/sum(!is.na(x))})
      for (i in 1:length(G_SNV_MAF)){
        if (G_SNV_MAF[i] <= 0.5){
          G_SNV[is.na(G_SNV[,i]),i] <- 0
        }
        else{
          G_SNV[is.na(G_SNV[,i]),i] <- 2
        }
      }

      pos_adj <- as.integer(seqGetData(genofile, "position"))
      ref_adj <- as.character(seqGetData(genofile, "$ref"))
      alt_adj <- as.character(seqGetData(genofile, "$alt"))
      variant_adj_info <- data.frame(chr,pos_adj,ref_adj,alt_adj)
      colnames(variant_adj_info) <- c("chr","pos","ref","alt")
      variant_adj_info

      seqResetFilter(genofile)
    }
  }

  summary_stat <- NULL
  cov_cond <- list()

  for(j in 1:ceiling((end_loc - start_loc + 1) / subsegment.size)){
    region_start_loc <- start_loc + (j-1) * subsegment.size
    region_end_loc <- start_loc + j * subsegment.size - 1
    region_end_loc <- min(region_end_loc,end_loc)
    print(paste0("Start Location: ", region_start_loc, ", End Location: ", region_end_loc, " (Subsegment: ", j, ")"))

    is.in <- (SNVlist)&(position>=region_start_loc)&(position<=region_end_loc)
    seqSetFilter(genofile,variant.id=variant.id[is.in],sample.id=phenotype.id)

    pos <- as.integer(seqGetData(genofile, "position"))
    ref <- as.character(seqGetData(genofile, "$ref"))
    alt <- as.character(seqGetData(genofile, "$alt"))
    if(check_qc_label){
      qc_label <- as.character(seqGetData(genofile, QC_label))
    }else{
      qc_label <- NULL
    }

    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))

    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index

    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]

    if(!is.null(Geno))
    {
      # Summary statistics
      genotype <- matrix_impute(Geno)
      variant_info <- data.frame(chr,pos,ref,alt,row.names=NULL)

      ## Annotation
      Anno.Int.PHRED.sub <- NULL
      Anno.Int.PHRED.sub.name <- NULL

      if(variant_type=="SNV")
      {
        if(Use_annotation_weights)
        {
          for(k in 1:length(Annotation_name))
          {
            if(Annotation_name[k]%in%Annotation_name_catalog$name)
            {
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
              Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

              if(Annotation_name[k]=="CADD")
              {
                Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
              }
              Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
            }
          }

          Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
          colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
        }
      }

      results_temp <- NULL
      try(results_temp <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                       Anno.Int.PHRED.sub,for_individual_analysis=TRUE),silent=silent)
      summary_stat <- rbind(summary_stat,results_temp)

      # Covariance matrices for conditional analysis
      if(!is.null(known_loci))
      {
        results_temp_list <- NULL
        try(results_temp_list <- MetaSTAARlite_worker_cov_cond(genotype,G_SNV,obj_nullmodel,variant_info,variant_adj_info),silent=silent)
        cov_cond$GTPG_cond <- rbind(cov_cond$GTPG_cond,results_temp_list$GTPG_cond)
        cov_cond$G_condTPG_cond <- results_temp_list$G_condTPG_cond
        cov_cond$variant_info <- rbind(cov_cond$variant_info,results_temp_list$variant_info)
        cov_cond$variant_adj_info <- results_temp_list$variant_adj_info
      }
    }
  }

  seqResetFilter(genofile)

  if(!is.null(known_loci))
  {
    return(list(summary_stat=summary_stat,
                cov_cond=cov_cond))
  }else
  {
    return(summary_stat)
  }
}
