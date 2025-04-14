#' Generates summary statistics of noncoding functional categories using MetaSTAARlite
#'
#' This function uses MetaSTAARlite to generate variant-level summary statistics
#' and sparse covariance matrices for noncoding functional categories of a gene.
#' @param chr an integer which specifies the chromosome number.
#' @param gene_name a character which specifies the name of the gene to be meta-analyzed using MetaSTAARlite.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function
#' in the \code{\link{STAARpipeline}} package, or the output from \code{fitNullModel} function in the \code{GENESIS} package
#' and transformed using the \code{\link{genesis2staar_nullmodel}} function in the \code{\link{STAARpipeline}} package.
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
#' @return a list of the following objects corresponding to each noncoding functional category of the given gene:
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
noncoding_MetaSTAARlite_worker <- function(chr,gene_name,genofile,obj_nullmodel,known_loci=NULL,
                                           cov_maf_cutoff=0.05,signif.digits=NULL,
                                           QC_label="annotation/filter",check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                           Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                           Use_annotation_weights=TRUE,Annotation_name=NULL,
                                           silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)

  phenotype.id <- as.character(obj_nullmodel$id_include)
  samplesize <- length(phenotype.id)

  #####################################
  #   Gene Info
  ## get SNV id
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
      
      ## variant id
      variant.id_adj <- seqGetData(genofile, "variant.id")
      ## get AF, Missing rate
      AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
      REF_AF <- AF_AC_Missing$af
      Missing_rate <- AF_AC_Missing$miss
      seqResetFilter(genofile)
      
      Genotype_sp <- Genotype_flip_sp_extraction(genofile,variant.id=variant.id_adj,
                                                 sample.id=phenotype.id,
                                                 REF_AF=REF_AF,
                                                 Missing_rate=Missing_rate,
                                                 QC_label=QC_label)
      G_SNV <- Genotype_sp$Geno
      
      if(!is.null(G_SNV) & inherits(G_SNV, "dgCMatrix"))
      {
        # geno_missing_imputation: "minor"
        G_SNV <- na.replace.sp(G_SNV,is_NA_to_Zero=TRUE)
        results_information_adj <- Genotype_sp$results_information
        variant_adj_info <- results_information_adj[,c("CHR","position","REF","ALT")]
        colnames(variant_adj_info) <- c("chr","pos","ref","alt")
        variant_adj_info
      }
      rm(Genotype_sp,AF_AC_Missing,REF_AF,Missing_rate,results_information_adj,variant.id_adj)
      gc()
    }
  }

  rm(filter)
  gc()

  ########################################
  #   Downstream

  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  is.in <- (GENCODE.Category=="downstream")&(SNVlist)
  variant.id.downstream <- variant.id[is.in]

  seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)

  rm(variant.id.downstream)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

  rm(GENCODE.Info)
  gc()

  rm(variant_gene_num)
  gc()

  Gene <- as.character(unlist(GENCODE.Info.split))

  rm(GENCODE.Info.split)
  gc()

  seqResetFilter(genofile)

  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

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

  summary_stat_list <- list()
  GTSinvG_rare_list <- list()
  cov_cond_list <- list()

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

  summary_stat_list[["downstream"]] <- summary_stat
  GTSinvG_rare_list[["downstream"]] <- GTSinvG_rare
  cov_cond_list[["downstream"]] <- cov_cond

  ########################################
  #   Upstream

  is.in <- (GENCODE.Category=="upstream")&(SNVlist)
  variant.id.upstream <- variant.id[is.in]

  seqSetFilter(genofile,variant.id=variant.id.upstream,sample.id=phenotype.id)

  rm(variant.id.upstream)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

  rm(GENCODE.Info)
  gc()

  rm(variant_gene_num)
  gc()

  Gene <- as.character(unlist(GENCODE.Info.split))

  rm(GENCODE.Info.split)
  gc()

  seqResetFilter(genofile)

  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

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

  summary_stat_list[["upstream"]] <- summary_stat
  GTSinvG_rare_list[["upstream"]] <- GTSinvG_rare
  cov_cond_list[["upstream"]] <- cov_cond

  ########################################################
  #                UTR

  is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
  variant.id.UTR <- variant.id[is.in]

  rm(GENCODE.Category)
  gc()

  seqSetFilter(genofile,variant.id=variant.id.UTR,sample.id=phenotype.id)

  rm(variant.id.UTR)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")

  rm(GENCODE.Info)
  gc()

  Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))

  rm(GENCODE.Info.split)
  gc()

  variant.id.SNV <- seqGetData(genofile, "variant.id")

  seqResetFilter(genofile)

  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

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

  summary_stat_list[["UTR"]] <- summary_stat
  GTSinvG_rare_list[["UTR"]] <- GTSinvG_rare
  cov_cond_list[["UTR"]] <- cov_cond

  #############################################
  #   Promoter-CAGE

  ## Promoter
  varid <- seqGetData(genofile, "variant.id")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

  # Subsetting Promoters that within +/-3kb of TSS and have CAGE signals
  CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
  CAGEBvt <- CAGEAnno!=""
  CAGEidx <- which(CAGEBvt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[CAGEidx])
  seqSetFilter(genofile,promGobj,intersect=TRUE)
  CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  CAGEGene <- unlist(lapply(strsplit(CAGEpromgene,"\\(|\\,|;|-"),`[[`,1))
  ##obtain variants info
  CAGEvchr <- as.numeric(seqGetData(genofile,"chromosome"))
  CAGEvpos <- as.numeric(seqGetData(genofile,"position"))
  CAGEvref <- as.character(seqGetData(genofile,"$ref"))
  CAGEvalt <- as.character(seqGetData(genofile,"$alt"))
  dfPromCAGEVarGene <- data.frame(CAGEvchr,CAGEvpos,CAGEvref,CAGEvalt,CAGEGene)

  ## get SNV id
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
  variant.id.SNV <- variant.id[SNVlist]

  dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
  dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
  dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
  dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)

  seqResetFilter(genofile)

  rm(dfPromCAGEVarGene)
  gc()

  ### Gene
  is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

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

  summary_stat_list[["promoter_CAGE"]] <- summary_stat
  GTSinvG_rare_list[["promoter_CAGE"]] <- GTSinvG_rare
  cov_cond_list[["promoter_CAGE"]] <- cov_cond

  ##################################################
  #       Promoter-DHS

  # Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
  rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
  rOCRsBvt <- rOCRsAnno!=""
  rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[rOCRsidx])

  seqSetFilter(genofile,promGobj,intersect=TRUE)
  rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
  ## obtain variants info
  rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
  rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
  rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
  rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
  dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)

  ## get SNV id
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
  variant.id.SNV <- variant.id[SNVlist]

  dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
  dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
  dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
  dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)

  seqResetFilter(genofile)

  rm(dfPromrOCRsVarGene)
  gc()

  ### Gene
  is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

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

  summary_stat_list[["promoter_DHS"]] <- summary_stat
  GTSinvG_rare_list[["promoter_DHS"]] <- GTSinvG_rare
  cov_cond_list[["promoter_DHS"]] <- cov_cond

  ###########################################
  #        Enhancer-CAGE

  #Now extract the GeneHancer with CAGE Signal Overlay
  genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  genehancer <- genehancerAnno!=""

  CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
  CAGE <- CAGEAnno!=""
  CAGEGeneHancervt <- CAGEAnno!=""&genehancerAnno!=""
  CAGEGeneHanceridx <- which(CAGEGeneHancervt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[CAGEGeneHanceridx])

  # variants that covered by whole GeneHancer without CAGE overlap.
  genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
  enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
  enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
  enhancervpos <- as.numeric(seqGetData(genofile,"position"))
  enhancervref <- as.character(seqGetData(genofile,"$ref"))
  enhancervalt <- as.character(seqGetData(genofile,"$alt"))
  dfHancerCAGEVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

  ## get SNV id
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
  variant.id.SNV <- variant.id[SNVlist]

  dfHancerCAGEVarGene.SNV <- dfHancerCAGEVarGene[SNVlist,]
  dfHancerCAGEVarGene.SNV$enhancervpos <- as.character(dfHancerCAGEVarGene.SNV$enhancervpos)
  dfHancerCAGEVarGene.SNV$enhancervref <- as.character(dfHancerCAGEVarGene.SNV$enhancervref)
  dfHancerCAGEVarGene.SNV$enhancervalt <- as.character(dfHancerCAGEVarGene.SNV$enhancervalt)

  seqResetFilter(genofile)

  rm(dfHancerCAGEVarGene)
  gc()

  ### Gene
  is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

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

  summary_stat_list[["enhancer_CAGE"]] <- summary_stat
  GTSinvG_rare_list[["enhancer_CAGE"]] <- GTSinvG_rare
  cov_cond_list[["enhancer_CAGE"]] <- cov_cond

  ##################################################
  #       Enhancer-DHS

  rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
  rOCRs <- rOCRsAnno!=""
  rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
  rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])
  # variants that covered by whole GeneHancer without rOCRs overlap.

  genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
  enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
  enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
  enhancervpos <- as.numeric(seqGetData(genofile,"position"))
  enhancervref <- as.character(seqGetData(genofile,"$ref"))
  enhancervalt <- as.character(seqGetData(genofile,"$alt"))
  dfHancerrOCRsVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

  rm(varid)
  gc()

  ## get SNV id
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
  variant.id.SNV <- variant.id[SNVlist]

  dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
  dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
  dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
  dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)

  seqResetFilter(genofile)

  rm(dfHancerrOCRsVarGene)
  gc()

  ### Gene
  is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

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

  summary_stat_list[["enhancer_DHS"]] <- summary_stat
  GTSinvG_rare_list[["enhancer_DHS"]] <- GTSinvG_rare
  cov_cond_list[["enhancer_DHS"]] <- cov_cond

  if(!is.null(known_loci))
  {
    return(list(summary_stat_list=summary_stat_list,
                GTSinvG_rare_list=GTSinvG_rare_list,
                cov_cond_list=cov_cond_list))
  }else
  {
    return(list(summary_stat_list=summary_stat_list,
                GTSinvG_rare_list=GTSinvG_rare_list))
  }
}
