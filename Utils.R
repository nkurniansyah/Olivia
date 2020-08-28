
readConfig <- function(file, ...) {
  config.table <- read.table(file, as.is=TRUE, ...)
  if (any(duplicated(config.table[, 1]))) stop("duplicated parameters in config file are not allowed!")
  config <- config.table[,2]
  names(config) <- config.table[,1]
  # recode tabs
  config[config %in% "\\t"] <- "\t"
  
  return(config)
}



writeConfig <- function(config, file, ...) {
  write.table(config, file=file, col.names=FALSE, ...)
}

#' @param required character vector of required parameter names
#' @param optional named vector of optional parameter values
#' @rdname readConfig
#'
#' @export
#'
setConfigDefaults <- function(config, required, optional) {
  # optional is a named list of default values
  default <- unname(optional)
  optional <- names(optional)
  
  config.params <- names(config)
  found.params <- intersect(config.params, c(required, optional))
  if (length(found.params) > 0) {
    message("found parameters: ", paste(found.params, collapse=", "))
  }
  
  # if required params not in config, stop
  missing.params <- setdiff(required, config.params)
  if (length(missing.params) > 0) {
    stop("missing required parameters: ", paste(missing.params, collapse=", "))
  }
  
  # if not in config, set default value
  set.params <- setdiff(optional, config.params)
  if (length(set.params) > 0) {
    config[set.params] <- default[match(set.params, optional)]
    message("using default values: ", paste(set.params, collapse=", "))
  }
  
  # note unsed params in config
  extra.params <- setdiff(config.params, c(required, optional))
  if (length(extra.params) > 0) {
    message("unused parameters: ", paste(extra.params, collapse=", "))
  }
  
  # return config with default values set
  config <- config[c(required, optional)]
  return(config)
}




#' Title: Emperical Pvalues
#'
#' @param p_values 
#' @param null_p_values 
#'
#' @return vector of emperical pvalues

quantile_emPval <- function(p_values, null_p_values){
  
  Fn <- ecdf(null_p_values)
  emp_pvalues <- Fn(p_values)
  
  inds <- which(emp_pvalues == 0)
  if (length(inds) > 0){
    emp_pvalues[inds] <- 1/(2*length(null_p_values))
  }
  
  emp_pvalues
}


#' Title: Normalization
#'
#' @param count_matrix : Row count (After basic filter)
#'
#' @return countsx: Normalize Row caunt
#' @export
#' 
normalized<- function(count_matrix){
  countsx <- t(t(count_matrix)/(colSums(count_matrix))*median(colSums(count_matrix)))
  return(countsx)
}


#' Title: Residual permutation for continous trait
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
#' @return value gene power with residual pemutations

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


#' Title: Log Transform method (Replcae zero wit minimum value divided by 2 to avoid inf for log 0)
#'
#' @param count_matrix : raw counts after gene filter 
#'
#' @return : log transfor gene counts


log_replace_half_min<- function(count_matrix){
  imputed_mat <- t(apply(count_matrix,1,function(x){x[x==0] <- min(x[x>0]/2);x}))
  imputed_mat <- log2(imputed_mat)
  return(imputed_mat)
}


#' Title:  Log Transform method (Add 0.5 to entire matrix)
#'
#' @param count_matrix  raw counts after gene filter 
#'
#' @return : log transfor gene counts
#' 

log_add_0.5<- function(count_matrix){
  imputed_mat <- log2(count_matrix + 0.5)
  return(imputed_mat)
}



#' Title: Log Transform method (Add min value divided by 2 to entire matrix)
#'
#' @param count_matrix  raw counts after gene filter 
#'
#' @return : log transfor gene counts
#' 
log_add_min <-  function(count_matrix){
  imputed_mat = t(apply(count_matrix,1,function(x){x= x + min(x[x>0]/2);x}))
  imputed_mat <- log2(imputed_mat)
  return(imputed_mat)
}

