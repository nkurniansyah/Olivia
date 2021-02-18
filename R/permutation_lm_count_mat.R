#' Fast linear regression
#'
#' Calculated associtation single gene count with multiple residual permutation as outcome
#'
#' @param residual_permutation Matrix of residual of permuation
#' @param covariates_string A character string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age,as.factor(sex)"
#' @param single_transcript  A vector of single transcript
#' @param pheno Phenotype data, includes the trait and covariates.
#' @return Linear regression results of  single gene permutation a data frame with columns permute(number of permute), beta,se,t_stat (t-statistic),p_value (permutation)
#' @examples
#' data(phenotype)
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' data("ENSG00000000003")
#' covariates_string<-"Age+Sex"
#' trait<-"Trait.1"
#' n_permute<- 100
#' permuted_trait<- sapply(seq_len(n_permute), function(x){
#'               permute_resids_trait(pheno = phenotype,
#'               trait = trait,
#'               covariates_string = covariates_string)})
#' lm_count_mat_permute(residual_permutation=permuted_trait ,covariates_string=covariates_string,
#'                       pheno=phenotype ,single_transcript= ENSG00000000003)
#' @export



lm_count_mat_permute<-function(residual_permutation, covariates_string, pheno, single_transcript){


  Ys<-as.matrix(residual_permutation)

  covariates_string<- as.character(covariates_string)

  #covars<- unlist(strsplit(covariates_string, ","))

  #model_string<- paste(covars, collapse = "+")


  model_string <- paste(covariates_string, "+", trait)
  XX<-model.matrix(as.formula(paste0("~", model_string)), data=pheno)
  stopifnot(nrow(XX) == nrow(Ys))
  stopifnot(nrow(XX) == length(single_transcript))

  # add single_transcript, and an intercept, to the design matrix
  XX <- cbind(single_transcript, XX)
  head(XX)

  numExplan <-ncol(XX)

  XtXinv <- solve(t(XX) %*% as.matrix(XX))
  XtXinv_se_arg <- sqrt(XtXinv["single_transcript","single_transcript"])

  XXproj <- XtXinv %*% t(XX)
  betas_mat <- XXproj %*% Ys
  betas <- betas_mat["single_transcript",]

  resid_Ys <-Ys - XX %*% XXproj %*% Ys
  sum_squares_resids <- colSums(resid_Ys^2)
  sigmas_square <- sum_squares_resids/(nrow(Ys)-numExplan)
  se_betas <- sqrt(sigmas_square)*XtXinv_se_arg

  test_stats <- betas/se_betas

  t_pval <- 2*pt(abs(test_stats), lower.tail=FALSE, df = nrow(Ys) - numExplan)
  res <- data.frame(permute = paste0("perm_",1:dim(Ys)[2]), beta = betas,
         se_beta = se_betas, test_stat = test_stats, pvalue = t_pval)

  return(res)
}




#' Wrapper function for permutation differential expression analysis
#'
#' The function is to calculate DEG (Differential Expression Genes) analysis for selected genes using multiple residual permuation results as outcome to calucate permutation p-value
#'
#' @param count_matrix A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#' @param pheno A data frame of phenotype data, includes the trait and covariates.
#' @param trait A character, the name of the exposure variable. The trait should be a column in pheno.
#' @param covariates_string A character string with specifying the covariats, include "as.factor" statements. example: covariate_string = "Age+as.factor(Sex)"
#' @param log_transform One of the transformations log_replace_half_min, log_add_min, log_add_0.5, or NULL (default)
#' @param gene_IDs A vector of selection of geneID, NULL if all genes are tested
#' @param n_permute Number of permutation. Default is 100000 times
#' @param seed Random seed
#' @param outcome_type continous and binary. Default is continous
#' @return Linear regression results as a data frame with columns GeneID, beta,se,t_stat (t-statistic), t_stat_df(degree of freedom),p_value, perm_pval
#' @examples
#' library(dplyr)
#' data(rnaseq_count_matrix)
#' data(phenotype)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' genes<-c("ENSG00000000003","ENSG00000000005","ENSG00000000419","ENSG00000000457","ENSG00000000460")
#' covariates_string<-"Age+Sex"
#' trait<-"Trait.1"
#' lm_count_mat_perm_pval(count_matrix=rnaseq_count_matrix, pheno=phenotype,
#'                        trait=trait, covariates_string= covariates_string,
#'                        gene_IDs=genes,n_permute=1000,
#'                        log_transform = "log_replace_half_min",
#'                        seed = NULL)
#' @export
#'



lm_count_mat_perm_pval <-function(count_matrix, pheno, trait, covariates_string,
                                 n_permute=100000,
                                 gene_IDs=NULL,
                                 log_transform = "log_replace_half_min",
                                 seed = NULL,
                                 outcome_type="continous"){

  if(is.null(gene_IDs)) message("No list gene ID/s are found. It will run permutations for all the genes and it will take long times")

  if (!is.null(seed)) set.seed(seed)

  deg <- lm_count_mat(count_matrix=count_matrix,
                      pheno=pheno,
                      trait=trait,
                      covariates_string=covariates_string,
                      gene_IDs=gene_IDs,
                      log_transform = log_transform)
  
  deg<- deg %>% dplyr::select(-fdr_bh)


  # generated permuted traits (by residual permutation)
  message("Performing residual permutation to generate permuted trait...")
  permuted_trait<- sapply(seq_len(n_permute), function(x){
  permute_resids_trait(pheno = pheno,
                      trait = trait,
                      covariates_string = covariates_string, outcome_type = outcome_type)
  })

  count_matrix <- log_transform_count(count_matrix, transform = log_transform)

  permute_pval<-lapply(seq_along(gene_IDs), function(y){
    transcript_val<- count_matrix[rownames(count_matrix)==gene_IDs[y],]
    transcript_df <- data.frame(transcript=transcript_val, stringsAsFactors = FALSE)
    rownames(transcript_df)<-names(transcript_val)
    current_pheno <- merge(pheno, transcript_df, by="row.names")
    head(current_pheno)
    permute_transcript <- lm_count_mat_permute(residual_permutation=permuted_trait,
                                               covariates_string=covariates_string,
                                               pheno=current_pheno,
                                               single_transcript=transcript_val)


    permute_pval<-permute_transcript$pvalue
    stopifnot(length(permute_pval)==nrow(transcript_val))



    deg_selected<-deg[which(deg$geneID==gene_IDs[y]),]
    perm_pval <- permutation_pvalues(pvalue = deg_selected$p_value, null_pval = permute_pval)
    new_deg<-cbind(deg_selected,perm_pval)

  })

  complete_deg<- do.call(rbind,permute_pval)
  return(complete_deg)
}
