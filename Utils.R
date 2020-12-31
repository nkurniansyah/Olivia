

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
#' @param phenotype :  Data frame of phenotype
#' @param outcome :  Character outcome, example: "ahi"
#' @param covariates_string :  Character covariates to adjust into model, example : "age,bmi,sex"

#' @return SizeFactor_normalization : Matrix of gene expression counts after normalization
#' 
#' 

SizeFactor<- function(count_matrix,covariates_string, outcome,phenotype){
  covariates_string<- as.character(covariates_string)
  designs<- gsub(",","+",covariates_string)
  
  des_matrix<- DESeqDataSetFromMatrix(countData = count_matrix, 
                                      colData = phenotype,
                                      design = formula(paste0("~ ",designs,"+",outcome)))
  
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

TMM<- function(count_matrix){
  
  counts <- DGEList(count_matrix)
  
  # normlize data using TMM method
  dgList<- calcNormFactors(counts, method = "TMM")
  
  TMM_normalization<- cpm(dgList)
  return(TMM_normalization)
}



#' Title: Log Transform method (Replcae zero wit minimum value divided by 2)
#'
#' @param count_matrix : A p x n matrix of gene expression counts. p are genes, n are individuals. Rownames are gene names
#'
#' @return : matrix of gene expression counts after transformation


log_replace_half_min<- function(count_matrix){
  imputed_mat <- t(apply(count_matrix,1,function(x){x[x==0] <- min(x[x>0]/2);x}))
  imputed_mat <- log2(imputed_mat)
  return(imputed_mat)
}


#' Title:  Log Transform method (Add 0.5 to entire matrix)
#'
#' @param count_matrix  A p x n matrix of gene expression counts. p are genes, n are individuals. Rownames are gene names
#'
#' @return : matrix of gene expression counts after transformation
#' 

log_add_0.5<- function(count_matrix){
  imputed_mat <- log2(count_matrix + 0.5)
  return(imputed_mat)
}



#' Title: Log Transform method (Add min value divided by 2 to entire matrix gene expression counts)
#'
#' @param count_matrix  A p x n matrix of gene expression counts. p are genes, n are individuals. Rownames are gene names
#'
#' @return : matrix of gene expression counts after transformation
#' 
log_add_min <-  function(count_matrix){
  imputed_mat = t(apply(count_matrix,1,function(x){x= x + min(x[x>0]/2);x}))
  imputed_mat <- log2(imputed_mat)
  return(imputed_mat)
}



#' Title: Emperical Pvalues
#'
#' @param p_values 
#' @param null_p_values 
#'
#' @return vector of emperical pvalues

quantile_empirical_pvalue <- function(p_values, null_p_values){
  
  Fn <- ecdf(null_p_values)
  emp_pvalues <- Fn(p_values)
  
  inds <- which(emp_pvalues == 0)
  if (length(inds) > 0){
    emp_pvalues[inds] <- 1/(2*length(null_p_values))
  }
  
  emp_pvalues
}



storey_empirical_pvalue <- function(z_score, null_z_score){
  emp_pvalues<- empPvals(z.score, null_zscore_res)
  emp_pvalues
}



#' Title: Residual permutation for continuous trait
#'
#' @param pheno : phenogtype
#' @param outcome : trait (example : "avgsat5, minsat5, rdi3p5)
#' @param covariates_string :cov to adjust 
#' @param seed : number of permutation
#'
#' @return (values with residual permuted)


permute_resids_trait <- function(pheno, outcome, covariates_string, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  # run the linear model
  fit <- lm(as.formula(paste(outcome, "~", covariates_string)), data = pheno)
  # get residual 
  resids <- fit$resid
  #create sample permutaions
  resid_permuted <- sample(resids)
  #get resiual + fit values
  resid_permuted + fitted.values(fit)	
}



#' Title: Residual permutaion for Binary trait
#'
#' @param pheno : phenogtype
#' @param outcome example : OSA
#' @param covariates_string :cov to adjust 
#' @param seed: number of permutation
#' @param family : "binomial"
#'
#' @return (values with residual permuted)

permute_resids_trait_bino <- function(pheno, outcome, covariates_string, seed = NULL, family){
  if (!is.null(seed)) set.seed(seed)
  # run the linear model
  fit<- glm(as.formula(paste(outcome , "~", covariates_string)), family = family,data = pheno)
  #fit <- lm(as.formula(paste(outcome, "~", covariates_string)), data = pheno)
  # get residual 
  rbinom(length(fit$fitted.values), size =1, prob = fit$fitted.values)
  #resids <- fit$resid
  #create sample permutaions
  #resid_permuted <- sample(resids)
  #get resiual + fit values
  #resid_permuted + fitted.values(fit)	
}


#rank-normalization function
rankNorm <- function(x){
  qnorm((rank(x)-0.5)/length(x))
}




#' Title: Power Simulations
#'
#' @param pheno : Phenotype
#' @param outcome : outcome : trait ("avgsat5, minsat5, rdi3p5)
#' @param covariates_string :cov to adjust 
#' @param gene_exp : selected gene to run the power
#' @param required_cor : corralations values (0.3,0.4,0.5)
#' @param seed : number of permutation
#'
#' @return value gene power with residual permutations

permute_resids_trait_cor <- function(pheno, outcome, covariates_string, gene_exp, required_cor, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  # run the linear model
  fit <- lm(as.formula(paste(outcome, "~", covariates_string)), data = pheno)
  # get residual 
  resids <- fit$resid
  # rank them in order to match with gene expression
  rank_gene_exp <- rank(gene_exp , ties.method = "random")
  rank_resid <- rank(resids,  ties.method = "random")
  # match them
  perfectly_matched_resids <- resids[order(resids)][rank_gene_exp]
  # create sample for permutations
  inds_to_scramble <- sample(1:length(rank_resid), size = floor((1-required_cor)*length(rank_resid)))
  resid_permuated <- perfectly_matched_resids
  # re-match resid permuation with index inds_to_scramble
  resid_permuated[inds_to_scramble] <- sample(resid_permuated[inds_to_scramble] )
  # return the resid permutaion and fit values 
  resid_permuated + fitted.values(fit)	
}
