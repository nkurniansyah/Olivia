

#' Title: Residual permutation for continuous trait
#'
#' @param pheno : a data.frame of phenotypes/variables
#' @param trait : a continuous variable to permute while preserving association with covariates
#' @param covariates_string :covariates to adjust for, as a string used in regression model
#' @param seed : random seed
#'
#' @return trait values after residual permuted


permute_resids_trait <- function(pheno, trait, covariates_string, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  # run the linear model
  fit <- lm(as.formula(paste(trait, "~", covariates_string)), data = pheno)
  # get residual 
  resids <- fit$resid
  #create sample permutations
  resid_permuted <- sample(resids)
  #get residual + fit values
  resid_permuted + fitted.values(fit)	
}



#' Title: Residual permutation for Binary trait
#'
#' @param pheno a data.frame of phenotypes/variables
#' @param trait a binary variable to be permuted while preserving assoication with covariates
#' @param covariates_string :covariates to adjust for, as a string provided to regression function
#' @param seed random seed
#' @param family : "binomial"
#'
#' @return permuted trait

permute_resids_trait_bino <- function(pheno, trait, covariates_string, seed = NULL, family){
  if (!is.null(seed)) set.seed(seed)
  # run the linear model
  fit <- glm(as.formula(paste(trait , "~", covariates_string)), family = family,data = pheno)
  # compute random trait with the same probabilities as the 
  # true traits as a function of covariates only
  rbinom(length(fit$fitted.values), size =1, prob = fit$fitted.values)
}


#rank-normalization function
rankNorm <- function(x){
  qnorm((rank(x)-0.5)/length(x))
}




#' Title: Power Simulations
#'
#' @param pheno : a data.frame with phenotypes
#' @param trait : a continuous variable to permute while preserving association with covariates
#' @param covariates_string :covariates to adjust for, as a string used in regression model
#' @param gene_exp : selected gene/transcript to simulate as associated with the trait
#' @param required_cor : corelations value simulating the strength of association between trait and gene_exp
#' @param seed : random seed
#'
#' @return permuted trait that is associated with gene_exp

permute_resids_trait_cor <- function(pheno, 
                                     trait, 
                                     covariates_string, 
                                     gene_exp, 
                                     required_cor, 
                                     seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  # run the linear model
  fit <- lm(as.formula(paste(trait, "~", covariates_string)), data = pheno)
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
  # re-match resid permutation with index inds_to_scramble
  resid_permuated[inds_to_scramble] <- sample(resid_permuated[inds_to_scramble] )
  # return the resid permutation and fit values 
  resid_permuated + fitted.values(fit)	
}
