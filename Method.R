
#' Title Linear Regressions
#'
#' @param count_matrix : Raw count (Transcript)
#' @param covariates_string : cova to adjust : format :  "Shipment,gender1,race1c,age5c,bmi5c,JHU,UMN"----- > NO SPACE
#' @param trait : outcome
#' @param pheno : phenotype data
#' @param log_transform  : log tarnsform : log_replace_min, log_add_min, log_add_0.5
#' @param list_geneID : vector of selection of geneID
#' @return Linear regrassion results
#' @export
#'
#' @examples
#' 

run_linear_regression <-function(count_matrix, pheno, trait, covariates_string, normal_pval = F, list_geneID=NULL, log_transform){
  
  #test_var<- trait
  if(!is.null(list_geneID)){
    count_matrix<-count_matrix[rownames(count_matrix) %in% list_geneID,]
    dim(count_matrix)}
  
  log_trans_mat<- log_transform(count_matrix)
  log_count_matrix<-t(log_trans_mat)
  
  stopifnot(nrow(log_count_matrix) == nrow(pheno))
  
  #stopifnot(rownames(log_count_matrix)pheno$ID))
  
  
  trait<- as.character(trait)
  #trait<- paste0("as.numeric(",trait,")")
  covars<- gsub(","," ",covariates_string)
  covars<- unlist(strsplit(covars,split = " "))
  
  #covariates_string <- paste0("as.factor(GENDER)+ as.integer(AGE)+as.factor(CENTER)+as.factor(gengrp6)+ ", trait)
  #print(covariates_string)
  #model1<- svyglm(as.formula(paste0(outcome, "~",trait,"+",paste(covars,collapse= "+"))),design=survey.design)
  
  cov <-c(covars, trait)
  # 
  XX<-model.matrix(as.formula(paste0("~", cov,collapse= "+")), data=pheno)
  XtXinv <- solve(t(XX) %*% as.matrix(XX))
  XtXinv_se_arg <- sqrt(XtXinv[trait,trait])
  numExplan <-ncol(XX)
  
  XXproj <- XtXinv %*% t(XX)
  betas_mat <- XXproj %*% log_count_matrix
  betas <- betas_mat[trait,]
  
  resid_Ys <-log_count_matrix - XX %*% XXproj %*% log_count_matrix
  sum_squares_resids <- colSums(resid_Ys^2)
  sigmas_square <- sum_squares_resids/(nrow(log_count_matrix)-numExplan)
  se_betas <- sqrt(sigmas_square)*XtXinv_se_arg
  
  test_stats <- betas/se_betas
  if (normal_pval){
    z_pval <- pchisq(test_stats^2, df = 1, lower.tail = FALSE)
    res <- data.frame(outcome = colnames(log_count_matrix), beta = betas, se_beta = se_betas, test_stat = test_stats, pvalue = z_pval)
  } else{ # p-value by t distribution
    t_pval <- 2*pt(abs(test_stats), lower.tail=FALSE, df = nrow(log_count_matrix) - numExplan)
    res <- data.frame(outcome = colnames(log_count_matrix), beta = betas, se_beta = se_betas, test_stat = test_stats, pvalue = t_pval)
  }
  head(res)
  
  colnames(res)<- c("GeneID", "Beta","SE","Stat","Pvalue")
  rownames(res)<-NULL 
  res<- res %>% mutate(FDR_BH= p.adjust(Pvalue, method = "BH"), z.score= qnorm(1-(Pvalue/2))*sign(Beta))
  return(res)
}



#' Title: Run Multivariates linear regression
#'
#' @param count_matrix : Raw count (Transcript)
#' @param pheno  phenotype data
#' @param covariates_string 
#' @param exposures : multiple exposusure to test : "avgsat5,minsat5,rdi3p5" ---> NO SPACE
#' @param normal_pval : P (based on Tstat)
#' @param list_geneID : vector of selection of geneID
#'
#' @return linear resgressions results
#' @export
#'
#' @examples
#' 
run_multivar_Linear_reggressions<- function(count_matrix, pheno, covariates_string, exposures, list_geneID=NULL, log_transform){
  
  if(!is.null(list_geneID)){
    count_matrix<-count_matrix[rownames(count_matrix) %in% list_geneID,]
    dim(count_matrix)
  }else{
    count_matrix<- count_matrix
  }
  log_trans_mat<- log_transform(count_matrix)
  print(dim(log_trans_mat))
  log_count_matrix<-t(log_trans_mat)
  
  exposures<- as.character(exposures)
  #exposures<- strsplit(exposures,",")[[1]]
  dim(log_trans_mat)
  print(exposures)
  
  covars<- gsub(","," ",covariates_string)
  covars<- unlist(strsplit(covars,split = " "))
  
  
  ### from the point of having an "X" matrix 
  cov <-c(covars, exposures)
  # 
  XX<-model.matrix(as.formula(paste0("~", cov,collapse= "+")), data=pheno)
  #XX <- cbind(1, as.matrix(X))
  print("here")
  XtXinv <- solve(t(XX) %*% as.matrix(XX))
  XtXinv_var_arg <- solve(XtXinv[exposures,exposures])
  numExplan <-ncol(XX)
  
  XXproj <- XtXinv %*% t(XX)
  
  betas_mat <- XXproj %*% log_count_matrix
  
  # effect sizes: each column correspond to a different transcript
  betas <- betas_mat[exposures,]
  betas_val<- t(betas)
  colnames(betas_val)<- c(paste0("Beta:",exposures))
  
  resid_Ys <-log_count_matrix - XX %*% XXproj %*% log_count_matrix
  sum_squares_resids <- colSums(resid_Ys^2)
  sigmas_square <- sum_squares_resids/(nrow(log_count_matrix)-numExplan)
  Joint_stats_arg1 <- XtXinv_var_arg %*% betas
  Joint_stats_arg2 <- colSums(betas*Joint_stats_arg1 )
  Joint_stats <- Joint_stats_arg2/sigmas_square
  Joint.Pvalue <- pchisq(Joint_stats, df = length(exposures), lower.tail = FALSE)
  
  res<- data.frame(betas_val,Joint_stats,Joint.Pvalue) %>% mutate(Joint.FDR_BH= p.adjust(Joint.Pvalue, method = "BH")) %>% rownames_to_column(var="GeneID")
  return(res)
}
