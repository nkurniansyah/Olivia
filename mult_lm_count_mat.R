#' Title Fast linear regression of each of gene expression counts on multiple traits
#'
#' @param count_matrix : A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#' @param pheno : phenotype data, includes the trait and covariates.  
#' @param traits : A character, the name of the exposure variables. The traits should columns in pheno.
#' @param covariates_string : a character string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age + as.factor(sex)"
#' @param log_transform_type  : One of the transformations log_replace_min, log_add_min, log_add_0.5, or NULL (default)
#' @param gene_IDs : vector of selection of gene IDs, NULL if all genes are tested
#' @return Linear regression results as a data frame with columns Gene_ID, Beta,SE,T_stat (t-statistic),P_value, Z_score (transformed by p-value)
#' @export
#'
#' @examples
#' 
mult_lm_count_mat <- function(count_matrix, pheno, 
                              covariates_string, traits, gene_IDs=NULL, log_transform = NULL){
  
  stopifnot(colnames(count_matrix) == rownames(pheno), all(is.element(traits, colnames(pheno))))
  
  count_matrix <- as.matrix(count_matrix)
  
  # if gene_IDs are (or is) provided, filter count_matrix to the requested genes
  if(!is.null(gene_IDs)) {
    count_matrix <- filter_by_genes(count_matrix, gene_IDs)
  }
  count_matrix <- log_transform_count(count_matrix, transform = log_transform_type)
  # transpose the matrix of counts:
  count_matrix <- t(count_matrix)
  
  # model for generating model matrix
  model_string <- paste(covariates_string, "+", paste(traits, collapse  = "+"))
  

  XX<-model.matrix(as.formula(paste0("~", model_string)), data=pheno)
  XtXinv <- solve(t(XX) %*% as.matrix(XX))
  XtXinv_var_arg <- solve(XtXinv[traits,traits])
  numExplan <-ncol(XX)
  
  XXproj <- XtXinv %*% t(XX)
  
  betas_mat <- XXproj %*% log_count_matrix
  
  # effect sizes: each column correspond to a different transcript
  betas <- betas_mat[traits,]
  betas_val<- t(betas)
  colnames(betas_val)<- c(paste0("Beta:",traits))
  
  resid_Ys <-log_count_matrix - XX %*% XXproj %*% log_count_matrix
  sum_squares_resids <- colSums(resid_Ys^2)
  sigmas_square <- sum_squares_resids/(nrow(log_count_matrix)-numExplan)
  Joint_stats_arg1 <- XtXinv_var_arg %*% betas
  Joint_stats_arg2 <- colSums(betas*Joint_stats_arg1 )
  Joint_stats <- Joint_stats_arg2/sigmas_square
  Joint_p_value <- pchisq(Joint_stats, df = length(traits), lower.tail = FALSE)
  
  res<- data.frame(betas_vals,Joint_stat = Joint_stats,Joint_p_value = Joint_p_value) %>% 
          mutate(Joint_FDR_BH= p.adjust(Joint_p_value, method = "BH")) %>% 
          rownames_to_column(var="GeneID")
  return(res)
}
