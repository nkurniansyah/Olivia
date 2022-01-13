
#' Residual permutation
#'
#' @param pheno A data.frame of phenotypes/variables
#' @param trait A continuous/binary variable to permute while preserving association with covariates
#' @param covariates_string Covariates to adjust for, as a string used in regression model
#' @param seed  Random seed
#' @param outcome_type continuous and binary, default is continuous
#' @return Vector of trait values after residual permuted
#' @examples
#' data(phenotype)
#' trait<-"Trait.1"; covars="Age+Sex"
#' permute_resids_trait(pheno=phenotype, trait=trait,
#'                      covariates_string=covars, 
#'                      seed = NULL, outcome_type="continuous")
#' @export



permute_resids_trait <- function(pheno, trait, covariates_string, seed = NULL, outcome_type="continuous"){

  if (!is.null(seed)) set.seed(seed)
  
  if(!all(trait %in% colnames(pheno))) stop("trait not found in the phenotype")


  stopifnot(trait %in% colnames(pheno))

  if (!is.element(outcome_type, c("continuous", "binary"))){
    stop(paste("Requested family type is", outcome_type,
               "allowed values are continuous and binary", "\n"))
  }

  

  if(outcome_type=="continuous"){
    fit<- lm(as.formula(paste(trait, paste(covariates_string),sep = "~")), data = pheno)

    # return fitted value + permuted residuals
    return(fitted.values(fit) + sample(fit$resid))
    
  }else if(outcome_type=="binary"){

    stopifnot(all(na.omit(pheno[,trait]) %in% 0:1))

    pheno[,trait]<-as.factor( pheno[,trait])
    fit<- glm(as.formula(paste(trait, paste(covariates_string),sep = "~")), family = binomial(link="logit"),data = pheno)
    
    # obtain estimated outcome probabilities
    estimated_probs <- fitted.values(fit)
    
    # return sampled outcomes according to the estimated probabilities
    return(rbinom(n = length(estimated_probs), size = 1, prob = estimated_probs))
  }
 
}


#' Power Simulations
#'
#' @param pheno A data.frame with phenotypes
#' @param trait A continuous/ binary variable to permute while preserving association with covariates
#' @param covariates_string Covariates to adjust for, as a string used in regression model
#' @param gene_exp Selected gene/transcript to simulate as associated with the trait
#' @param required_cor Corelations value simulating the strength of association between trait and gene_exp
#' @param seed Random seed
#' @param outcome_type continuous and binary
#' @return A vector permuted trait that is associated with gene_exp
#' @examples
#' data(phenotype)
#' data(ENSG00000000003)
#' trait<-"Trait.1"; covars="Age+Sex"
#' permute_resids_trait_cor(pheno= phenotype,trait=trait,
#'                          covariates_string=covars, required_cor=0.3,
#'                          gene_exp=ENSG00000000003,seed=NULL, outcome_type="continuous")
#' @export



permute_resids_trait_cor <- function(pheno,
                                     trait,
                                     covariates_string,
                                     gene_exp,
                                     required_cor,
                                     seed = NULL,
                                     outcome_type="continuous"){

  if (!is.null(seed)) set.seed(seed)



  if(!all(trait %in% colnames(pheno))) stop("trait not found in the phenotype")


  stopifnot(trait %in% colnames(pheno))

  if (!is.element(outcome_type, c("continuous", "binary"))){
    stop(paste("Requested family type is", outcome_type,
               "allowed values are continuous and binary", "\n"))
  }

  if(outcome_type=="continuous"){

    fit<- lm(as.formula(paste(trait, paste(covariates_string),sep = "~")), data = pheno)

  }else if(outcome_type=="binary"){

    stopifnot(all(na.omit(pheno[,trait]) %in% 0:1))
    pheno[,trait]<-as.factor( pheno[,trait])
    fit<- glm(as.formula(paste(trait, paste(covars),sep = "~")), family = binomial(link="logit"),data = pheno)

  }

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
