#' Merges the generated summary statistics and sparse covariance matrices of different studies
#'
#' This function merges the generated variant-level summary statistics and sparse covariance matrices
#' of all participating studies in preparation for the meta-analysis step of MetaSTAARlite.
#' @param chr an integer which specifies the chromosome number.
#' @param sample.sizes a numeric vector with the length of \code{study.names}
#' indicating the sample size of each study.
#' @param sumstat.list a list containing study-specific summary statistics from all participating studies.
#' @param cov.list a list containing study-specific sparse weighted covariance matrices from all participating studies.
#' @param rare_maf_cutoff a numeric value specifying the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param cov_maf_cutoff a numeric vector with the length of \code{study.names}
#' indicating the maximum minor allele frequency cutoffs under which the sparse weighted
#' covariance files between variants are stored.
#' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label} (default = TRUE).
#' @param variant_type a character value specifying the type of variant included in the analysis. Choices include
#'  "SNV", "Indel", or "variant" (default = "SNV").
#' @param Use_annotation_weights a logical value which determines if annotations will be used as weights or not (default = TRUE).
#' @param Annotation_name a character vector of annotation names used in MetaSTAARlite (default = NULL).
#' @return a list with the following members:
#' @return \code{info}: the merged data frame of all variants in the variant-set
#' of interest whose combined minor allele frequency is below \code{rare_maf_cutoff}, including the
#' following information (listed in the same order as \code{U} and the rows/columns of \code{cov}):
#' chromosome (chr), position (pos), reference allele (ref), alternative allele (alt),
#' combined minor allele count (MAC), and combined minor allele frequency (MAF).
#' @return \code{U} the merged score statistics vector of all variants in the variant-set
#' of interest whose combined minor allele frequency is below \code{rare_maf_cutoff}.
#' @return \code{cov} the merged covariance matrix of all variants in the variant-set
#' of interest whose combined minor allele frequency is below \code{rare_maf_cutoff}.
#' @return \code{annotation_phred} the merged functional annotation data in PHRED score scale
#' of all variants in the variant-set of interest whose combined minor allele frequency is below \code{rare_maf_cutoff}.
#' @references Li, X., et al. (2023). Powerful, scalable and resource-efficient
#' meta-analysis of rare variant associations in large whole genome sequencing studies.
#' \emph{Nature Genetics}, \emph{55}(1), 154-164.
#' (\href{https://doi.org/10.1038/s41588-022-01225-6}{pub})
#' @export

MetaSTAARlite_merge <- function(chr,sample.sizes,sumstat.list,cov.list,
                                rare_maf_cutoff=0.01,cov_maf_cutoff,
                                check_qc_label=TRUE,variant_type=c("SNV","Indel","variant"),
                                Use_annotation_weights=TRUE,Annotation_name=NULL){

  ## evaluate choices
  variant_type <- match.arg(variant_type)

  cov_maf_cutoff[cov_maf_cutoff == 0.5] <- 0.5 + 1e-16
  ### summary statistics
  sumstat.list <- lapply(sumstat.list, function(x) {
    if (!(is.null(x)) && !("qc_label" %in% colnames(x))){
      data.frame(x[,1:4],qc_label="PASS",x[,5:dim(x)[2]],stringsAsFactors = FALSE)
    }else{
      x
    }
  })
  sumstat.varid.list <- mapply(function(x,y) {
    x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),1:4]
  }, x = sumstat.list, y = cov_maf_cutoff, SIMPLIFY = FALSE)
  sumstat.varid.merge <- do.call("rbind",sumstat.varid.list)
  sumstat.varid.nodup <- sumstat.varid.merge[!duplicated(sumstat.varid.merge),]
  if (is.null(sumstat.varid.nodup)) {
    return(NULL)
  }else if (dim(sumstat.varid.nodup)[1] == 0) {
    return(NULL)
  }
  sumstat.merge.list <- lapply(sumstat.list, function(x) {
    if (is.null(x)) {
      cbind(sumstat.varid.nodup,qc_label=NA,alt_AC=NA,MAC=NA,MAF=NA,N=NA,U=NA,V=NA,"1"=NA)
    }else {
      left_join(sumstat.varid.nodup,x,by=c("chr"="chr",
                                           "pos"="pos",
                                           "ref"="ref",
                                           "alt"="alt"))
    }
  })
  sumstat.merge.list <- mapply(function(x,y) {
    x[is.na(x[,"N"]),"N"] <- y
    x[is.na(x[,"qc_label"]),"qc_label"] <- "PASS"
    x[is.na(x)] <- 0
    return(x)
  }, x = sumstat.merge.list, y = sample.sizes, SIMPLIFY = FALSE)
  sumstat.varid.nodup$index <- 1:dim(sumstat.varid.nodup)[1]
  sumstat.index.list <- lapply(sumstat.varid.list, function(x) {
    if (is.null(x)) {
      integer(0)
    }else{
      left_join(x,sumstat.varid.nodup,by=c("chr"="chr",
                                           "pos"="pos",
                                           "ref"="ref",
                                           "alt"="alt"))$index
    }
  })
  rm(sumstat.list,sumstat.varid.list,sumstat.varid.merge)
  gc()

  ### covariance matrices
  cov.list <- lapply(cov.list, function(x) {
    if (is.null(x)) {
      as(matrix(nrow=0,ncol=0),"dgCMatrix")
    }else if (dim(x)[1] == 0) {
      as(matrix(nrow=0,ncol=0),"dgCMatrix")
    }else {
      x[,1:dim(x)[1]]
    }
  })

  cov.null <- matrix(0, nrow = dim(sumstat.varid.nodup)[1], ncol = dim(sumstat.varid.nodup)[1])
  cov.merge.list <- mapply(function(x,y) {
    if (length(x) == 1) {
      # in this case y is a scalar
      cov.null[x,x] <- y
    }else if (length(x) > 1) {
      cov.null[x,x] <- as.matrix(forceSymmetric(y, uplo = "U"))
    }
    return(cov.null)
  }, x = sumstat.index.list, y = cov.list, SIMPLIFY = FALSE)

  ### select rare variant based on the input cutoff
  alt_AC.merge <- as.integer(Reduce("+",lapply(sumstat.merge.list, function(x) {x$alt_AC})))
  N.merge.nonzero <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N * (x$MAC != 0)}))
  alt_AF.merge.nonzero <- alt_AC.merge / (2 * N.merge.nonzero)
  N.merge.zero <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N * (x$MAC == 0)}))
  alt_AC.merge <- alt_AC.merge + (alt_AF.merge.nonzero > 0.5) * (2 * N.merge.zero)
  N.merge <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N}))
  MAC.merge <- pmin(alt_AC.merge, 2 * N.merge - alt_AC.merge)
  MAF.merge <- MAC.merge / (2 * N.merge)
  rv.index <- (MAF.merge<rare_maf_cutoff) & Reduce("*",mapply(function(x,y) {
    x$MAF<y
  }, x = sumstat.merge.list, y = cov_maf_cutoff, SIMPLIFY = FALSE))
  if (check_qc_label){
    rv.index <- rv.index & Reduce("*",lapply(sumstat.merge.list, function(x) {
      x$qc_label=="PASS"
    }))
  }

  info <- cbind(sumstat.varid.nodup[,c("chr","pos","ref","alt")],
                MAC=MAC.merge,MAF=MAF.merge)[rv.index,]

  U.merge <- Reduce("+", lapply(sumstat.merge.list, function(x) {
    x <- x[rv.index,]
    return(x$U)
  }))

  cov.merge <- Reduce("+", mapply(function(x,y) {
    (x - as.matrix(y[,(which(colnames(y)=="V")+1):dim(y)[2]]) %*% t(as.matrix(y[,(which(colnames(y)=="V")+1):dim(y)[2]])))[rv.index, rv.index]
  }, x = cov.merge.list, y = sumstat.merge.list, SIMPLIFY = FALSE))

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
        Annotation.list <- lapply(sumstat.merge.list, function(x) {
          if(Annotation_name[k]%in%colnames(x)) {
            x[[Annotation_name[k]]]
          }else {
            0
          }
        })
        Annotation.PHRED.sum <- Reduce("+", Annotation.list)
        Annotation.PHRED.gt0 <- Reduce("+",lapply(Annotation.list, function(x) {x > 0}))
        Annotation.PHRED <- Annotation.PHRED.sum / Annotation.PHRED.gt0
        Annotation.PHRED[!is.finite(Annotation.PHRED)] <- 0

        if(Annotation_name[k]=="aPC.LocalDiversity")
        {
          Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
          Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
        }
        Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- c(Anno.Int.PHRED.sub.name)
    }
  }

  annotation_phred <- Anno.Int.PHRED.sub[rv.index,]

  rm(list=setdiff(ls(), c("info","U.merge","cov.merge","annotation_phred")))
  gc()

  return(list(info=info,
              U=U.merge,
              cov=as.matrix(cov.merge),
              annotation_phred=annotation_phred))
}
