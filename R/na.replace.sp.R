#' Replace Missing Values in a Sparse Genotype Matrix
#'
#' The \code{na.replace.sp} function replaces missing values (NA) in a sparse genotype matrix 
#' (`dgCMatrix` format). If `is_NA_to_Zero = TRUE`, NA values are replaced with 0. 
#' Otherwise, NA values in each column are replaced with the corresponding entries in `m`.
#' This function is inspired by \code{glmnet::na.replace} for sparse matrices.
#'
#' @param genotype_sp A sparse genotype matrix of class `dgCMatrix` from the \code{Matrix} package.
#' @param m A numeric vector specifying the replacement values for each column.
#' @param is_NA_to_Zero A logical value indicating whether NA values should be replaced with 0 
#' (default: FALSE). If FALSE, NA values are replaced column-wise using `m`.
#' 
#' @return A `dgCMatrix` object with missing values replaced accordingly.
#' 
#' @examples
#' library(Matrix)
#' set.seed(123)
#' # Create a sparse matrix with some NA values
#' mat <- Matrix(c(1, NA, 3, 0, NA, 2, 4, 5, NA), nrow = 3, sparse = TRUE)
#' print(mat)
#' 
#' # Replace NA values with 0
#' mat_imputed <- na.replace.sp(mat, m = c(0.5, 1, 1.5), is_NA_to_Zero = TRUE)
#' print(mat_imputed)
#'
#' # Replace NA values with values from m
#' mat_imputed_m <- na.replace.sp(mat, m = c(0.5, 1, 1.5), is_NA_to_Zero = FALSE)
#' print(mat_imputed_m)
#'
#' @export

na.replace.sp <- function(genotype_sp, m, is_NA_to_Zero = FALSE) 
{
  if (!inherits(genotype_sp, "dgCMatrix")) stop("genotype_sp must be a dgCMatrix")
  
  x <- genotype_sp@x
  i <- genotype_sp@i
  p <- genotype_sp@p
  ncol <- ncol(genotype_sp)
  
  nas <- is.na(x)
  
  if (any(nas)) 
  {
    if (is_NA_to_Zero) 
    {
      # Remove NA values (set NA to 0)
      x_new <- x[!nas]
      i_new <- i[!nas]
      
      # Compute new column pointers (p array)
      nnz_per_col_new <- sapply(seq_len(ncol), function(j) 
      {
        if (p[j] < p[j + 1]) 
        {
          sum(!nas[(p[j] + 1):p[j + 1]])  # Count non-NA values per column
        } else 
        {
          0
        }
      })
      
      p_new <- as.integer(c(0, cumsum(nnz_per_col_new)))
      
      # Update the sparse matrix
      genotype_sp@x <- x_new
      genotype_sp@i <- i_new
      genotype_sp@p <- p_new
      
    } else 
    {
      # Replace NA values with m[j], where j is the column index
      ccount <- diff(p)  # Number of nonzero elements per column
      cindex <- rep(seq_len(ncol), ccount)
      genotype_sp@x[nas] <- m[cindex[nas]]
    }
  }
  return(genotype_sp)
}
