#' Title Fast linear regression of each of gene expression counts on an exposure
#'
#' @param count_matrix : A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#' @param pheno : phenotype data, includes the trait and covariates.  
#' @param trait : A character, the name of the exposure variable. The trait should be a column in pheno.
#' @param covariates_string : a character string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age,as.factor(sex)"
#' @param log_transform_type  : One of the transformations log_replace_min, log_add_min, log_add_0.5, or NULL (default)
#' @param gene_IDs : vector of selection of geneID, NULL if all genes are tested
#' @return Linear regression results as a data frame with columns GeneID, Beta,SE,T_stat (t-statistic),P_value, Z_score (transformed by p-value)
#' @export
#'
#' @examples
#' 

lm_count_mat <-function(count_matrix, pheno, trait, covariates_string,
                                  gene_IDs=NULL, log_transform = NULL){
  
  stopifnot(colnames(count_matrix) == rownames(pheno$ID), is.element(trait, colnames(pheno)))
  
  count_matrix <- as.matrix(count_matrix)
  
  # if gene_IDs are (or is) provided, filter count_matrix to the requested genes
  if(!is.null(gene_IDs)) {
    count_matrix <- filter_by_genes(count_matrix, gene_IDs)
  }
  
  count_matrix <- log_transform_count(count_matrix, transform = log_transform)

  # transpose the matrix of counts:
  count_matrix <- t(count_matrix)
  
  covariates_string<- as.character(covariates_string)
  cov<- gsub(",","+",covariates_string)
  model_string <- paste(cov, "+", trait)
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
  res <- data.frame(geneID = colnames(count_matrix), beta = betas, 
                      se = se_betas, t_stat = test_stats, p_value = t_pval)
  
  rownames(res)<-NULL 
  res<- res %>% mutate(FDR_BH= p.adjust(p_value, method = "BH"), 
                        z_score= qnorm(1-(p_value/2))*sign(beta))
  return(res)
}








#' Title Wrapper function for differential expression analysis, includes computation of empirical p-values
#'
#' @param count_matrix : A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#' @param pheno : phenotype data, includes the trait and covariates.  
#' @param trait : A character, the name of the exposure variable. The trait should be a column in pheno.
#' @param covariates_string : a character string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age,as.factor(sex)"
#' @param log_transform_type  : One of the transformations log_replace_min, log_add_min, log_add_0.5, or NULL (default)
#' @param gene_IDs : vector of selection of geneID, NULL if all genes are tested
#' @return Linear regression results as a data frame with columns GeneID, Beta,SE,T_stat (t-statistic),P_value, Z_score (transformed by p-value)
#' @export
#'
#' @examples
#' 

lm_count_mat_emp_pval <-function(count_matrix, pheno, trait, covariates_string, 
                                 n_permute=100,
                                 gene_IDs=NULL, 
                                 log_transform = NULL, 
                                 seed = NULL,
                                 empirical_type = "storey"){
  
  if (!is.null(seed)) set.seed(seed)
  
  deg <- lm_count_mat(count_matrix=count_matrix, 
                      pheno=pheno, 
                      trait=trait, 
                      covariates_string=covariates_string,
                      gene_IDs=gene_IDs, 
                      log_transform = log_transform)

  permuted_trait <- matrix(NA, nrow = nrow(pheno), ncol = n_permute)
  
  # generated permuted traits (by residual permutation)
  message("Performing residual permutation to generate permuted trait...")
  for (i in 1:n_permute){
    permuted_trait[,i] <- permute_resids_trait(pheno = phenotypes,
                                               trait = trait,
                                               covariates_string = covariates_string)
  }
  
  # perform differential expression analysis for each of the permuted traits and save z-scores
  null_z_scores <-  matrix(NA, nrow = nrow(deg), ncol = n_permute)
  
  message(paste("performing differential expression analysis on", n_permute, "permuted traits"))
  for (i in 1:n_permute){
    pheno$perm_trait <- permuted_trait[,i]
    deg_perm <- lm_count_mat(count_matrix=count_matrix, 
                             pheno=pheno, 
                             trait="perm_trait", 
                             covariates_string=covariates_string,
                             gene_IDs=gene_IDs, 
                             log_transform = log_transform)
    null_z_scores[,i] <- deg_perm$z_score
  }
  
  message("Computing empirical p-values")
  emp_pvals <- compute_empirical_pvalues(statistics = deg$z_score,
                            null_statistics = as.numeric(null_z_scores),
                            stat_type = "z_score", 
                            empirical_type = empirical_type)
  res <- deg
  res$emp_pvals <- emp_pvals
  res$BH_emp_pvals <- p.adjust(res$emp_pvals, method = "BH")
  return(res)
}
