

#' Title Normalize a matrix of transcript counts
#'
#' @param count_matrix : A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#' @param normalized_type  : One of the normalization method median_normalization, SizeFactor, TMM, or NULL 
#' @return  A p x n matrix of gene expression counts after normalization
#' @export
#' 

normalize_trancript_count <- function(count_matrix, normalized_type = "median_normalization"){
  if (is.null(median_normalization)){
    message("No normalization method of gene counts requested")
    return(count_matrix)
  }
  if (!is.element(normalized_type, c("median_normalization", "SizeFactor", "TMM"))){
    stop("Requested normalization not allowed. 
         Allowed normalization names are median_normalization, SizeFactor, TMM")
  }
  
  if (normalized_type == "median_normalization") {
    message("Applying median normalization  ...")
    
    return(median_normalization(count_matrix))
  }
  
  if (normalized_type == "SizeFactor"){
    message("Applying size factor ...")
    
    return(SizeFactor(count_matrix))
  }
  
  if (normalized_type == "TMM"){
    message("Applying TMM ...")
    
    return(TMM(count_matrix))
  }
  
  
  }

if(normalize_type=="median_libsize"){
  message("applying median library size ...")
  count_matrix <- normalize_median(count_matrix) 
  
}else if (normalize_type=="SizeFactor"){
  message("Applying EstimeSizeFactor ...")
  #normalize_sizefactor<- function(count_matrix, covariates_string, outcome, phenotype){
  
  count_matrix <- normalize_sizefactor(count_matrix, covariates_string = covariates_string, outcome =trait, phen=pheno) 
  
}else if(normalize_type=="TMM"){
  message("Applying TMM  ...")
  
  
  count_matrix <- normalize_tmm(count_matrix) 
  
}


