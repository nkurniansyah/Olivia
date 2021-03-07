### residual plot


#' Title
#'
#' @param traits Trait/s, the name of the exposure variable. The trait/s should be a column in pheno.
#' @param pheno A data frame of phenotype, includes the trait and covariates
#' @param covariates_string Characters string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age,as.factor(sex)"
#' @return residual plot
#'
#' @examples
#' library(ggplot2)
#' library(reshape2)
#' data(phenotype)
#' covariates_string<-"Age+Sex"
#' traits<-c("Trait.1","Trait.2")
#' residual_plot(pheno=phenotype, covariates_string = covariates_string,traits=traits)
#' @export



residual_plot<- function(pheno, covariates_string, traits){
  
  stopifnot(is.element(traits, colnames(pheno)))

  message(paste(c("Generate residual of ",traits,"...")), collapse= "  ")
  resid_out<-list()
  
  for(trait in traits){
    fit<-lm(as.formula(paste(trait, paste(covariates_string), sep = " ~ ")), data=pheno)
    
    #fit<- lm(as.formula(paste(outcome , "~", cov_str)), data = com)
    resid<- data.frame(fit$residuals)
    colnames(resid)<-trait
    resid_out[[trait]]<-resid
  }
  
  all_resid<- do.call(cbind, resid_out)
  head(all_resid)
  all_resid$no<- rownames(all_resid)
  resid_traits<- melt(all_resid, id.vars = "no")
  head(resid_traits)
  
  message("Generate the residual plot...")
  resid_plot<-ggplot(resid_traits, aes(value, fill=variable ))+geom_density(alpha=0.3)+ theme_bw() +ggtitle(paste0("Residual Plots"))+guides(fill=FALSE) 
  
  return(resid_plot+ facet_wrap(~ variable, scales = "free"))
  

}







#' Title
#'
#' @param pheno  A data frame of phenotype, includes the trait and covariates
#' @param categorical_variable  all selected categorical variable in model example: c("Sex", "Race")
#' @param numeric_variable all selected numeric variable in model example: c("Age", "Trait.1", "Trait.2")
#' @strata strata Variable to stratified example: "Race"
#' @return sumary phenotype 
#' @export
#'
#' @examples
#' library(tableone)
#' data(phenotype)
#' races<- c("Asian","White","Black")
#' phenotype$Race<-sample(races, nrow(pheno), replace = T)

#' categorical_variable<-c("Sex","Race")
#' numeric_variable<-c("Age","Trait.1","Trait.2")

summary_phenotype(pheno=phenotype, categorical_variable =categorical_variable,  numeric_variable=numeric_variable, strata="Race")

#' 

summary_phenotype<- function(pheno, categorical_variable, numeric_variable, strata){
  
  stopifnot(is.element(categorical_variable, colnames(pheno)))
  
  stopifnot(is.element(numeric_variable, colnames(pheno)))
  

  all_variable<- c(categorical_variable,numeric_variable )
  
  message(paste(c("Generate summary phnotype using",all_variable ),collapse = " "))
  
  tableOne <- CreateTableOne(vars =all_variable ,data = pheno, factorVars = categorical_variable, test=F,strata =strata )
  
  return(print(tableOne))
  
}
