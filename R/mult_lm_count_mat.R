#' Fast linear regression for multiple exposure
#'
#' Calculated associtation gene counts matrix with the pehonotype where each of gene counts on multiple traits
#'
#' @param count_matrix  A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#' @param pheno A data frame of phenotype data, includes the trait and covariates.
#' @param traits Characters, the name of the exposure variables. The traits should columns in pheno.
#' @param covariates_string A character string with specifying the covariates, include "as.factor" statements. example: covariates_string = "age + as.factor(sex)"
#' @param log_transform One of the transformations log_replace_half_min, log_add_min, log_add_0.5, or NULL (default)
#' @param gene_IDs A vector of selection of gene IDs, NULL if all genes are tested
#' @return Linear regression results as a data frame with columns geneID, adjLogFC.Trait.1,adjLogFC.Trait.2 ,se,chisq_stat,chisq_stat_df,p_value(join p-value),fdr_bh
#' @examples
#' set.seed(123)
#' library(dplyr)
#' data(phenotype)
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' traits<-c("Trait.1","Trait.2")
#' covars<- "Age+Sex"
#' mult_lm_count_mat(count_matrix=rnaseq_count_matrix,pheno = phenotype,traits = traits,
#'                   covariates_string=covars)
#' @export
#'

mult_lm_count_mat <- function(count_matrix, pheno, covariates_string, traits,
                              gene_IDs=NULL, log_transform = "log_replace_half_min"){
  traits<- as.character(traits)


  stopifnot(colnames(count_matrix) == rownames(pheno), all(is.element(traits, colnames(pheno))))


  if(length(traits)<2) stop("Only found single trait, run linear regression for single trait instead")


  count_matrix <- as.matrix(count_matrix)
  # if gene_IDs are (or is) provided, filter count_matrix to the requested genes
  if(!is.null(gene_IDs)) {
    count_matrix <- filter_by_genes(count_matrix, gene_IDs)
  }
  count_matrix <- log_transform_count(count_matrix, transform = log_transform)
  # transpose the matrix of counts:
  count_matrix <- t(count_matrix)



  covariates_string<- as.character(covariates_string)



  model_string <- c(covariates_string, traits)
  model_string<- paste(model_string, collapse = "+")


  XX<-model.matrix(as.formula(paste0("~", model_string)), data=pheno)
  XtXinv <- solve(t(XX) %*% as.matrix(XX))
  XtXinv_var_arg <- solve(XtXinv[traits,traits])
  numExplan <-ncol(XX)

  XXproj <- XtXinv %*% t(XX)

  betas_mat <- XXproj %*% count_matrix

  # effect sizes: each column correspond to a different transcript
  betas <- betas_mat[traits,]
  betas_val<- t(betas)
  colnames(betas_val)<- c(paste0("adjLogFC_",traits))

  resid_Ys <-count_matrix - XX %*% XXproj %*% count_matrix
  sum_squares_resids <- colSums(resid_Ys^2)
  sigmas_square <- sum_squares_resids/(nrow(count_matrix)-numExplan)
  Joint_stats_arg1 <- XtXinv_var_arg %*% betas
  Joint_stats_arg2 <- colSums(betas*Joint_stats_arg1 )
  Joint_stats <- Joint_stats_arg2/sigmas_square
  Joint_p_value <- pchisq(Joint_stats, df = length(traits), lower.tail = FALSE)


  res<- data.frame(geneID = colnames(count_matrix),betas_val,chisq_stat = Joint_stats,
                   chisq_stat_df = length(traits), p_value = Joint_p_value) %>%
                   mutate(fdr_bh= p.adjust(p_value, method = "BH"))

  return(res)
}


#' Wrapper function for differential expression analysis for multiple exposure
#'
#' The function is to calculate DEG (Differential Expression Genes) analysis for multiple exposure using residual permuation approach to calculate empirical p-value
#'
#' @param count_matrix  A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#' @param pheno A data frame phenotype data, includes the trait and covariates.
#' @param traits Characters, the name of the exposure variable. The traits should be a column in pheno.
#' @param covariates_string A character string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age+as.factor(sex)"
#' @param log_transform One of the transformations log_replace_half_min, log_add_min, log_add_0.5, or NULL (default)
#' @param gene_IDs : vector of selection of geneID, NULL if all genes are tested
#' @param n_permute number of computing residual permutation. Default is 100 times
#' @param seed Random seed
#' @param outcome_type continuous and binary.Default is continuous
#' @return Linear regression results as a data frame with columns geneID, beta.Trait1,beta.Trait1 ,se,
#'        chisq_stat,chisq_stat_df(degree of freedom),p_value(join p-value),emp_pvals,bh_emp_pvals
#'
#' @examples
#' set.seed(123)
#' library(dplyr)
#' data(phenotype)
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' traits<-c("Trait.1","Trait.2")
#' covars<- "Age+Sex"
#' lm_mult_count_mat_emp_pval(count_matrix=rnaseq_count_matrix, pheno=phenotype, traits=traits, covariates_string=covars,
#'                         outcome_type="continuous")
#' @export
#'


lm_mult_count_mat_emp_pval <-function(count_matrix, pheno, traits, covariates_string,
                                      n_permute=100, gene_IDs=NULL,seed = NULL,
                                      log_transform = "log_replace_half_min",outcome_type="continuous"){

  if (!is.null(seed)) set.seed(seed)



  deg <- mult_lm_count_mat(count_matrix=count_matrix,
                      pheno=pheno,
                      traits=traits,
                      covariates_string=covariates_string,
                      gene_IDs=gene_IDs,
                      log_transform = log_transform)

  message("Performing residual permutation to generate permuted trait...")

  traits<- as.character(traits)

  multi_permuted_trait<- lapply(seq_along(traits),function(z){
    permuted_trait<-sapply(seq_len(n_permute), function(x){
      permute_resids_trait(pheno = pheno,
                           trait= traits[z], seed = seed,
                           covariates_string = covariates_string, outcome_type = outcome_type)
    })
  })

  names(multi_permuted_trait)<-traits

  message(paste("performing differential expression analysis on", n_permute, "permuted traits"))

  permute_val<-lapply(seq_len(n_permute), function(y){
    perm_trait<-sapply(multi_permuted_trait, `[`,, y)
    colnames(perm_trait)<- paste0("perm_",colnames(perm_trait))
    pheno<- cbind(pheno,perm_trait)
    head(pheno)
    multi_exp<- colnames(pheno)[grepl("perm_", colnames(pheno))]
    deg_perm <- mult_lm_count_mat(count_matrix=count_matrix,
                                  pheno=pheno,
                                  traits =multi_exp,
                                  covariates_string=covariates_string,
                                  gene_IDs=gene_IDs,
                                  log_transform = log_transform)
    null_stat<-deg_perm[,"p_value"]
  })


  null_pval<-unlist(permute_val)

  stopifnot(length(null_pval)==n_permute*nrow(deg))

  message("Computing quantile empirical p-values")

  emp_pvals <- compute_quantile_empirical_pvalues(statistics = deg[,"p_value"],
                                                  null_statistics = as.numeric(null_pval))
  deg$emp_pvals <- emp_pvals
  deg<- deg %>% mutate(bh_emp_pvals= p.adjust(emp_pvals, method = "BH"))
  rownames(deg)<-NULL
  
  return(deg)
}
