MetaSTAARlite_merge <- function(chr,sample.sizes,sumstat.list,cov.list,
                                rare_maf_cutoff=0.01,cov_maf_cutoff,
                                check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL){
  
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