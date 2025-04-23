#' Generates summary statistics of long noncoding RNA (ncRNA) category using MetaSTAARlite
#'
#' This function uses MetaSTAARlite to generate variant-level summary statistics
#' and sparse covariance matrices the exonic and splicing category of an ncRNA gene.
#' @param chr an integer which specifies the chromosome number.
#' @param gene_name a character which specifies the name of the gene to be meta-analyzed using MetaSTAARlite.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function
#' in the \code{STAARpipeline} package, or the output from \code{fitNullModel} function in the \code{GENESIS} package
#' and transformed using the \code{\link{genesis2staar_nullmodel}} function in the \code{STAARpipeline} package.
#' @param known_loci the data frame of variants to be adjusted for in conditional analysis. Should
#' contain four columns in the following order: chromosome (CHR), position (POS), reference allele (REF),
#' and alternative allele (ALT) (default = NULL).
#' @param cov_maf_cutoff a numeric value indicating the maximum minor allele frequency cutoff
#' under which the sparse weighted covariance file between variants is stored (default = 0.05).
#' @param signif.digits an integer indicating the number of significant digits to be used
#' for storing the sparse weighted covariance file. If \code{signif.digits} is NULL,
#' it is assumed that no rounding will be performed (default = NULL).
#' @param QC_label a character specifying the channel name of the QC label in the GDS/aGDS file
#' (default = "annotation/filter").
#' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label}.
#' If \code{check_qc_label} is FALSE, then the summary statistics will be stored for PASS variants from the study.
#' If \code{check_qc_label} is TRUE, then the summary statistics will be stored for all variants from the study,
#' together will an additional column of \code{qc_label} (default = FALSE).
#' @param variant_type a character value specifying the type of variant included in the analysis. Choices include
#'  "SNV", "Indel", or "variant" (default = "SNV").
#' @param Annotation_dir a character specifying the channel name of the annotations in the aGDS file
#' (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the annotation name and the corresponding
#' channel name in the aGDS file.
#' @param Use_annotation_weights a logical value which specifies if annotations will be used as weights
#' or not (default = TRUE).
#' @param Annotation_name a character vector of annotation names used in MetaSTAARlite (default = NULL).
#' @param silent a logical value which determines if the report of error messages will be suppressed (default = FALSE).
#' @return a list of the following objects corresponding to the exonic and splicing category of the given ncRNA gene:
#' (1) the data frame of all variants in the variant-set (the variant-level summary statistics file),
#' including the following information: chromosome (chr), position (pos), reference allele (ref),
#' alternative allele (alt), quality control status (qc_label, optional), alternative allele count (alt_AC), minor allele count (MAC),
#' minor allele frequency (MAF), study sample size (N), score statistic (U), variance (V), variant annotations specified in
#' \code{Annotation_name}, and the low-rank decomposed component of the covariance file;
#' (2) the sparse matrix of all variants in the variant-set whose minor allele frequency is below \code{cov_maf_cutoff} (the sparse weighted
#' covariance file); (3) the summary statistics and covariance matrices corresponding to the specified gene for variants to be conditioned on
#' in \code{known_loci}.
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
ncRNA_MetaSTAARlite_worker <- function(chr,gene_name,genofile,obj_nullmodel,known_loci=NULL,
                                       cov_maf_cutoff=0.05,signif.digits=NULL,
                                       QC_label="annotation/filter",check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                       Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                       Use_annotation_weights=TRUE,Annotation_name=NULL,
                                       silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)

  phenotype.id <- as.character(obj_nullmodel$id_include)
  samplesize <- length(phenotype.id)

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

  variant.id <- seqGetData(genofile, "variant.id")
  
  G_SNV <- NULL
  if(!is.null(known_loci))
  {
    position <- as.integer(seqGetData(genofile, "position"))
    allele <- seqGetData(genofile, "allele")

    loc_SNV <- c()
    for (i in 1:dim(known_loci)[1]){
      loc_SNV <- c(loc_SNV,which((position==known_loci$POS[i])&(allele==paste0(known_loci$REF[i],",",known_loci$ALT[i]))))
    }

    if(length(variant.id[loc_SNV])>=1)
    {
      seqSetFilter(genofile,variant.id=variant.id[loc_SNV],sample.id=phenotype.id)
      
      G_SNV <- seqGetData(genofile, "$dosage")
      
      if (!is.null(G_SNV)){
        
        ## genotype id
        id.genotype <- as.character(seqGetData(genofile,"sample.id"))
        
        id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
        phenotype.id.merge <- data.frame(phenotype.id)
        phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
        id.SNV.match <- phenotype.id.merge$index
        
        G_SNV <- G_SNV[id.SNV.match,,drop=FALSE]
        G_SNV <- matrix_impute(G_SNV)
        
        pos_adj <- as.integer(seqGetData(genofile, "position"))
        ref_adj <- as.character(seqGetData(genofile, "$ref"))
        alt_adj <- as.character(seqGetData(genofile, "$alt"))
        variant_adj_info <- data.frame(chr,pos_adj,ref_adj,alt_adj)
        colnames(variant_adj_info) <- c("chr","pos","ref","alt")
        variant_adj_info
      }
      seqResetFilter(genofile)
    }
  }

  rm(filter)
  gc()

  ## ncRNA SNVs
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  is.in <- ((GENCODE.Category=="ncRNA_exonic")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.Category=="ncRNA_splicing"))&(SNVlist)

  variant.id.ncRNA <- variant.id[is.in]

  rm(GENCODE.Category)
  gc()

  seqSetFilter(genofile,variant.id=variant.id.ncRNA,sample.id=phenotype.id)

  rm(variant.id.ncRNA)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
  Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))

  Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
  Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
  Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))

  rm(GENCODE.Info)
  gc()

  rm(GENCODE.Info.split)
  gc()

  variant.id.ncRNA <- seqGetData(genofile, "variant.id")

  seqResetFilter(genofile)

  ### Gene
  is.in <- union(which(Gene_list_1==gene_name),which(Gene_list_2==gene_name))
  is.in <- union(is.in,which(Gene_list_3==gene_name))

  variant.is.in <- variant.id.ncRNA[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## variant id
  variant.is.in <- seqGetData(genofile, "variant.id")
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(length(variant.is.in)>=1)
  {
    # Annotation
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
    
    ## get AF, Missing rate
    AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
    REF_AF <- AF_AC_Missing$af
    Missing_rate <- AF_AC_Missing$miss
    rm(AF_AC_Missing)
  } else
  {
    REF_AF <- Missing_rate <- NULL
  }
  
  Genotype_sp <- Genotype_flip_sp_extraction(genofile,variant.id=variant.is.in,
                                             sample.id=phenotype.id,
                                             REF_AF=REF_AF,
                                             Missing_rate=Missing_rate,
                                             annotation_phred=Anno.Int.PHRED.sub,
                                             QC_label=QC_label)
  Geno <- Genotype_sp$Geno
  Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
  results_information <- Genotype_sp$results_information
  rm(Genotype_sp)
  gc()
  
  if(check_qc_label){
    qc_label <- results_information$qc_label
  }else{
    qc_label <- NULL
  }

  summary_stat <- NULL
  GTSinvG_rare <- NULL
  cov_cond <- NULL

  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    variant_info <- results_information[,c("CHR","position","REF","ALT")]
    colnames(variant_info) <- c("chr","pos","ref","alt")
    
    ALT_AF <- results_information$ALT_AF
    
    # geno_missing_imputation: "minor"
    Geno <- na.replace.sp(Geno,is_NA_to_Zero=TRUE)
    
    ## Summary statistics
    try(summary_stat <- MetaSTAARlite_worker_sumstat(Geno,ALT_AF,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub),silent=silent)

    ## Covariance matrices
    try(GTSinvG_rare <- MetaSTAARlite_worker_cov(Geno,obj_nullmodel,cov_maf_cutoff,
                                                 qc_label,signif.digits),silent=silent)

    ## Covariance matrices for conditional analysis
    if(!is.null(known_loci) & !is.null(G_SNV))
    {
      try(cov_cond <- MetaSTAARlite_worker_cov_cond(Geno,G_SNV,obj_nullmodel,variant_info,variant_adj_info),silent=silent)
    }
  }

  seqResetFilter(genofile)

  if(!is.null(known_loci))
  {
    return(list(summary_stat=summary_stat,
                GTSinvG_rare=GTSinvG_rare,
                cov_cond=cov_cond))
  }else
  {
    return(list(summary_stat=summary_stat,
                GTSinvG_rare=GTSinvG_rare))
  }
}
