
#' Residual permutation
#'
#' @param pheno A data.frame of phenotypes/variables
#' @param trait A continuous/binary variable to permute while preserving association with covariates
#' @param covariates_string Covariates to adjust for, as a string used in regression model
#' @param seed  Random seed
#' @param family A description of the error distribution to be used in the model. gaussian if the variable is continous,
#'               and binomial if variable is binary. The default is "gaussian"
#' @return Vector of trait values after residual permuted
#' @examples
#' data(phenotype)
#' trait<-"Trait.1"; covars="Age,Sex"
#' permute_resids_trait(pheno=phenotype, trait=trait,
#'                      covariates_string=covars, seed = NULL, family="gaussian")
#' @export


permute_resids_trait <- function(pheno, trait, covariates_string, seed = NULL, family="gaussian"){

  if (!is.null(seed)) set.seed(seed)
  covars<- unlist(strsplit(covariates_string, ","))

  if(!all(covars %in% colnames(pheno))) stop("covariates not found in the phenotype")

  if(!all(trait %in% colnames(pheno))) stop("trait not found in the phenotype")


  stopifnot(trait %in% colnames(pheno))

  if (!is.element(family, c("gaussian", "binomial"))){
    stop(paste("Requested family type is", family,
               "allowed values are gaussian and binomial", "\n"))
  }

  #covariates_string<- as.character(covariates_string)
  #cov<- gsub(",","+",covariates_string)

  if(family=="gaussian"){
    fit<- lm(as.formula(paste(trait, paste(covars, collapse = "+"),sep = "~")), data = pheno)

  }else if(family=="binomial"){

    stopifnot(all(na.omit(pheno[,trait]) %in% 0:1))

    pheno[,trait]<-as.factor( pheno[,trait])
    fit<- glm(as.formula(paste(trait, paste(covars, collapse = "+"),sep = "~")), family = binomial(link="logit"),data = pheno)

  }
  # get residual
  resids <- fit$resid
  #create sample permutations
  resid_permuted <- sample(resids)
  #get residual + fit values
  resid_permuted + fitted.values(fit)
}


#' Power Simulations
#'
#' @param pheno A data.frame with phenotypes
#' @param trait A continuous/ binary variable to permute while preserving association with covariates
#' @param covariates_string Covariates to adjust for, as a string used in regression model
#' @param gene_exp Selected gene/transcript to simulate as associated with the trait
#' @param required_cor Corelations value simulating the strength of association between trait and gene_exp
#' @param seed Random seed
#' @param family A description of the error distribution to be used in the model. gaussian if the variable is continous,
#'               and binomial if variable is binary. The default is "gaussian"
#' @return A vector permuted trait that is associated with gene_exp
#' @examples
#' data(phenotype)
#' data(ENSG00000000003)
#' trait<-"Trait.1"; covars="Age,Sex"
#' permute_resids_trait_cor(pheno= phenotype,trait=trait,
#'                          covariates_string=covars, required_cor=0.3,
#'                          gene_exp=ENSG00000000003,seed=NULL, family="gaussian")
#' @export



permute_resids_trait_cor <- function(pheno,
                                     trait,
                                     covariates_string,
                                     gene_exp,
                                     required_cor,
                                     seed = NULL,
                                     family="gaussian"){

  if (!is.null(seed)) set.seed(seed)

  covars<- unlist(strsplit(covariates_string, ","))

  if(!all(covars %in% colnames(pheno))) stop("covariates not found in the phenotype")

  if(!all(trait %in% colnames(pheno))) stop("trait not found in the phenotype")


  stopifnot(trait %in% colnames(pheno))

  if (!is.element(family, c("gaussian", "binomial"))){
    stop(paste("Requested family type is", family,
               "allowed values are gaussian and binomial", "\n"))
  }


  if(family=="gaussian"){

    fit<- lm(as.formula(paste(trait, paste(covars, collapse = "+"),sep = "~")), data = pheno)

  }else if(family=="binomial"){

    stopifnot(all(na.omit(pheno[,trait]) %in% 0:1))
    pheno[,trait]<-as.factor( pheno[,trait])
    fit<- glm(as.formula(paste(trait, paste(covars, collapse = "+"),sep = "~")), family = binomial(link="logit"),data = pheno)

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
