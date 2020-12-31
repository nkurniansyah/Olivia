

#' Title Normalize a matrix of transcript counts
#'
#' @param count_matrix : A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#' @param normalized_type  : One of the normalization method median_normalization, SizeFactor, TMM, or NULL 
#' @return  A p x n matrix of gene expression counts after normalization
#' @export
#' 


normalize_trancript_count <- function(count_matrix,
                                      normalization_type = "median_normalization",
                                      phenotypes=NULL,
                                      covariates_string=NULL, 
                                      trait=NULL){
  if (is.null(normalization_type)){
    message("No normalization method of gene counts requested")
    return(count_matrix)
  }
  if (!is.element(normalization_type, c("median_normalization", "SizeFactor", "TMM"))){
    stop("Requested normalization not allowed. 
         Allowed normalization names are median_normalization, SizeFactor, TMM")
  }
  
  if (normalization_type == "median_normalization") {
    message("Applying median normalization...")
    
    return(median_normalization(count_matrix))
  }
  
  if (normalization_type == "SizeFactor"){
    message("Applying size factor normalization...")
    return(SizeFactor(count_matrix,
                      phenotypes = phenotypes,
                      covariates_string = covariates_string,
                      trait = trait))
  }
  
  if (normalization_type == "TMM"){
    message("Applying TMM normalization...")
    
    return(TMM(count_matrix))
  }
  
}



#' Title: Median Normalization
#'
#' @param count_matrix :  A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#'
#' @return median_normalization : Matrix of gene expression counts after normalization
#' 
#' 
median_normalization <- function(count_matrix){
  median_normalization <- t(t(count_matrix)/(colSums(count_matrix))*median(colSums(count_matrix)))
  return(median_normalization)
}



#' Title: SizeFactor Normalization
#'
#' @param count_matrix :  A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#' @param phenotypes :  Data frame of phenotype
#' @param trait :  Character trait of interest
#' @param covariates_string :  Character covariates to adjust into model, example : "age,bmi,sex"

#' @return SizeFactor_normalization : Matrix of gene expression counts after normalization
#' 
#' 

SizeFactor <- function(count_matrix, phenotypes, covariates_string, trait){
  
  des_matrix <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                      colData = phenotype,
                                      design = formula(paste0("~ ",covariates_string,"+",trait)))
  
  des_matrix <- estimateSizeFactors(des_matrix)
  
  SizeFactor_normalization <- counts(des_matrix, normalized=TRUE)
  
  return(SizeFactor_normalization)
}


#' Title: TMM Normalization
#'
#' @param count_matrix :  A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#'
#' @return TMM_normalization : Matrix of gene expression counts after normalization
#' 
#' 

TMM <- function(count_matrix){
  
  counts <- DGEList(count_matrix)
  
  # normlize data using TMM method
  dgList<- calcNormFactors(counts, method = "TMM")
  
  TMM_normalization<- cpm(dgList)
  return(TMM_normalization)
}


