#' Title Fast linear regression of each of gene expression counts on an exposure
#'
#' @param count_matrix : A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#' @param pheno : phenotype data, includes the trait and covariates.  
#' @param trait : A character, the name of the exposure variable. The trait should be a vector in pheno.
#' @param covariates_string : a character string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age + as.factor(sex)"
#' @param log_transform_type  : One of the transformations log_replace_min, log_add_min, log_add_0.5, or NULL (default)
#' @param geneIDs : vector of selection of geneID, NULL if all genes are tested
#' @return Linear regression results as a data frame with columns GeneID, Beta,SE,T_stat (t-statistic),P_value, Z_score (transformed by p-value)
#' @export
#'
#' @examples
#' 

lm_count_mat <-function(count_matrix, pheno, trait, covariates_string, normal_pval = FALSE, geneIDs=NULL, log_transform = NULL){
  
  stopifnot(colnames(count_matrix) == rownames(pheno), is.element(trait, colnames(pheno)))
  
  count_matrix <- as.matrix(count_matrix)
  
  # if geneIDs are (or is) provided, filter count_matrix to the requested genes
  # make sure that requested genes are in the data!
  if(!is.null(geneIDs)){
    if (!all(geneIDs %in% rownames(count_matrix))) {
      genes_not_avail <- geneIDs[which(!geneIDs %in% rownames(count_matrix))]
      message(paste("Some requested geneIDs not in count_matrix, proceeding without using ", 
                    paste(genes_not_avail, collapse = ", ")))
    }
    genes_avail <- geneIDs[which(geneIDs %in% rownames(count_matrix))]
    if (length(genes_avail) == 0){
      stop("None of the requested genes are in the count_matrix, stopping...")
    }
    message(paste0("Filtering count_matrix to genes", genes_avail))
    
    count_matrix<-count_matrix[rownames(count_matrix) %in% geneIDs, , drop = FALSE]
  }
  
  count_matrix <- log_transform_count(count_matrix, transform = log_transform_type)
  
  # transpose the matrix of counts:
  count_matrix <- t(count_matrix)
  
  model_string <- paste(covariates_string, "+", trait)
  # 
  XX<-model.matrix(as.formula(paste0("~", model_string)), data=pheno)
  XtXinv <- solve(t(XX) %*% as.matrix(XX))
  XtXinv_se_arg <- sqrt(XtXinv[trait,trait])
  numExplan <-ncol(XX)
  
  XXproj <- XtXinv %*% t(XX)
  betas_mat <- XXproj %*% count_matrix
  betas <- betas_mat[trait,]
  
  resid_Ys <- count_matrix - XX %*% XXproj %*% count_matrix
  sum_squares_resids <- colSums(resid_Ys^2)
  sigmas_square <- sum_squares_resids/(nrow(count_matrix)-numExplan)
  se_betas <- sqrt(sigmas_square)*XtXinv_se_arg
  
  test_stats <- betas/se_betas
  t_pval <- 2*pt(abs(test_stats), lower.tail=FALSE, df = nrow(count_matrix) - numExplan)
  res <- data.frame(outcome = colnames(count_matrix), beta = betas, se_beta = se_betas, test_stat = test_stats, pvalue = t_pval)
  
  colnames(res)<- c("GeneID", "Beta","SE","T_stat","P_value")
  rownames(res)<-NULL 
  res<- res %>% mutate(FDR_BH= p.adjust(Pvalue, method = "BH"), Z_score= qnorm(1-(P_value/2))*sign(Beta))
  return(res)
}
