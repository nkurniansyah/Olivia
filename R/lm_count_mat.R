#' Fast linear regression
#'
#' Calculated associtation gene counts matrix with the pehonotype where each of gene counts on aexposure
#'
#' @param count_matrix  A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#' @param pheno  A data frame of phenotype, includes the trait and covariates.
#' @param trait A character, the name of the exposure variable. The trait should be a column in pheno.
#' @param covariates_string Characters string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age,as.factor(sex)"
#' @param log_transform  One of the transformations log_replace_half_min, log_add_min, log_add_0.5, or NULL (default)
#' @param gene_IDs Vector of selection of geneID, NULL if all genes are tested
#' @return Linear regression results as a data frame with columns geneID, beta,se,t_stat (t-statistic),t_stat_df(degree of freedom),p_value,fdr_bh 
#' @examples
#' library(dplyr)
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' data(phenotype)
#' trait<-"Trait.1"
#' covars<-"Age+Sex"
#' log_transform<-"log_replace_half_min"
#' lm_count_mat(count_matrix=rnaseq_count_matrix,pheno=phenotype,trait=trait,
#'              covariates_string=covars, log_transform=log_transform)
#' @export
#'

lm_count_mat <-function(count_matrix, pheno, trait, covariates_string,
                        gene_IDs=NULL, log_transform = "log_replace_half_min"){

  stopifnot(colnames(count_matrix) == rownames(pheno), is.element(trait, colnames(pheno)))

  count_matrix <- as.matrix(count_matrix)

  # if gene_IDs are (or is) provided, filter count_matrix to the requested genes
  if(!is.null(gene_IDs)) {
    count_matrix <- filter_by_genes(count_matrix, gene_IDs)
  }

  count_matrix <- log_transform_count(count_matrix, transform = log_transform)

  count_matrix <- t(count_matrix)

  #covariates_string<- as.character(covariates_string)

  #covars<- unlist(strsplit(covariates_string, ","))

  model_string <- paste0(covariates_string,"+",trait)
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
  t_stat_df <- nrow(count_matrix) - numExplan
  t_pval <- 2*pt(abs(test_stats), lower.tail=FALSE, df = t_stat_df)
  res <- data.frame(geneID = colnames(count_matrix), beta = betas, se = se_betas, 
                    t_stat = test_stats,t_stat_df=t_stat_df, p_value = t_pval)

  rownames(res)<-NULL
  
  res<- res %>% mutate(fdr_bh= p.adjust(p_value, method = "BH"))
  return(res)
}


#' Wrapper function for differential expression analysis
#'
#' The function is to calculate DEG (Differential Expression Genes) analysis using residual permuation approach to calculate empirical p-value
#'
#' @param count_matrix A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#' @param pheno A data frame of phenotype data, includes the trait and covariates
#' @param trait A character, the name of the exposure variable. The trait should be a column in pheno
#' @param covariates_string Characters string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age+as.factor(sex)"
#' @param log_transform  One of the transformations log_replace_half_min, log_add_min, log_add_0.5, or NULL (default)
#' @param gene_IDs A vector of selection of geneID, NULL if all genes are tested
#' @param n_permute number of computing residual permutation. Default is 100 times
#' @param seed Random seed
#' @param outcome_type continous and binary.Default is continous
#' @return Linear regression results as a data frame with columns geneID, beta,se,t_stat (t-statistic),t_stat_df(degree of freedom),p_value, fdr_bh,,
#'         emp_pvals,bh_emp_pvals
#' @examples
#' set.seed(123)
#' library(dplyr)
#' data(phenotype)
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' trait<-"Trait.1"
#' covars<- "Age+Sex"
#' lm_count_mat_emp_pval(count_matrix=rnaseq_count_matrix,pheno = phenotype,trait = trait,
#'                       covariates_string=covars, outcome_type="continous")
#' @export


lm_count_mat_emp_pval <-function(count_matrix, pheno, trait, covariates_string,
                                 n_permute=100,
                                 gene_IDs=NULL,
                                 log_transform = "log_replace_half_min",
                                 seed = NULL,
                                 outcome_type="continous"){

  if (!is.null(seed)) set.seed(seed)

  deg <- lm_count_mat(count_matrix=count_matrix,
                      pheno=pheno,
                      trait=trait,
                      covariates_string=covariates_string,
                      gene_IDs=gene_IDs,
                      log_transform = log_transform)

  message("Performing residual permutation to generate permuted trait...")


  permuted_trait<-sapply(seq_len(n_permute), function(x){
                  permute_resids_trait(pheno = pheno,
                  trait= trait, seed = seed,
                  covariates_string = covariates_string,outcome_type = outcome_type)
  })

  message(paste("performing differential expression analysis on", n_permute, "permuted traits"))


  permute_val<-lapply(seq_len(n_permute), function(y){
    pheno$perm_trait <- permuted_trait[,y]
    deg_perm <- lm_count_mat(count_matrix=count_matrix,
                             pheno=pheno,
                             trait="perm_trait",
                             covariates_string=covariates_string,
                             gene_IDs=gene_IDs,
                             log_transform = log_transform)
    null_pval<-deg_perm[,"p_value"]
    
  })

  null_pval<-unlist(permute_val)
  head(null_pval)

  stopifnot(length(null_pval)==n_permute*nrow(deg))

  message("Computing quantile empirical p-values")

  emp_pvals <- compute_quantile_empirical_pvalues(statistics = deg[,"p_value"],
                                                  null_statistics = as.numeric(null_pval))
  
  deg$emp_pvals <- emp_pvals
  deg<- deg %>%mutate(bh_emp_pvals=p.adjust(emp_pvals, method = "BH"))
  
  return(deg)
}


