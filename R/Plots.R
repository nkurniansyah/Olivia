
#'
#' @param traits Trait/s, the name of the exposure variable(s). The trait/s should be a column in pheno.
#' @param pheno A data frame of phenotypes, includes the trait and covariates
#' @param covariates_string Characters string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age+as.factor(sex)"
#' @return residual plot
#'
#' @examples
#' library(ggplot2)
#' library(reshape2)
#' data(phenotype)
#' covariates_string<-"Age+Sex+Race"
#' traits<-c("Trait.1","Trait.2")
#' residual_plot(pheno=phenotype, covariates_string = covariates_string,traits=traits)
#' @export



residual_plot<- function(pheno, covariates_string, traits){
  
  stopifnot(is.element(traits, colnames(pheno)))

  message(paste(c("Generate residual of ",traits,"...")), collapse= "  ")
  resid_out<-list()
  
  for(trait in traits){
    fit<-lm(as.formula(paste(trait, paste(covariates_string), sep = " ~ ")), data=pheno)
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
  resid_plot<-ggplot(resid_traits, aes(value, fill=variable )) +
                  geom_density(alpha=0.3) + 
                  theme_bw() +
                  ggtitle(paste0("Residual Plots")) +
                  guides(fill=FALSE) 
  
  return(resid_plot+ facet_wrap(~ variable, scales = "free"))
  
}


#' Title Summarize phenotypes
#'
#' @param pheno  A data frame of phenotypes, includes the trait and covariates
#' @param categorical_variables  all selected categorical variables in the model, example: c("Sex", "Race")
#' @param numeric_variables all selected numeric variables in the model, example: c("Age", "Trait.1", "Trait.2")
#' @param strata Variable to stratify, to example: "Race"
#' @return summary table of phenotypes 
#' @export
#' @examples
#' library(tableone)
#' data(phenotype)
#' categorical_variable<-c("Sex","Race")
#' numeric_variable<-c("Age","Trait.1","Trait.2")
#' strata<-"Race"
#' summarize_phenotypes(pheno=phenotype, categorical_variables =categorical_variable,  numeric_variables=numeric_variable, strata=strata)

summarize_phenotypes <- function(pheno, categorical_variables, numeric_variables, strata){
  
  stopifnot(is.element(categorical_variables, colnames(pheno)))
  
  stopifnot(is.element(numeric_variables, colnames(pheno)))
  

  all_variables <- c(categorical_variables,numeric_variables )
  
  message(paste(c("Generate summary of phenotype using",all_variables ),collapse = " "))
  
  tableOne <- CreateTableOne(vars =all_variables,strata =strata ,data = pheno, factorVars = categorical_variables, test=F, )

  return(tableOne)
  
}




#' Title Violin plot
#'
#' @param pheno A data frame of phenotype, includes the trait and covariates
#' @param strata  Variable to stratified example: "Race"
#' @param norm_count_matrix A matrix of gene counts (possibly normalize transformed). rows are genes, columns are individuals
#' @param selected_transcript transcript selection example: "ENSG00000002549"
#' @param log_transform  One of the transformations log_replace_half_min, log_add_min, log_add_0.5, or NULL (default)
#'
#' @return violin plot 
#' @export
#'
#' @examples
#' strata<-"Race"
#' data(phenotype)
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' selected_transcript<- "ENSG00000002549"
#' strata<-"Race"
#' violin_plot(pheno=phenotype, strata=strata, norm_count_matrix=rnaseq_count_matrix, selected_transcript=selected_transcript )


violin_plot<- function(pheno, strata, norm_count_matrix, selected_transcript){
  
  stopifnot(is.element(strata, colnames(pheno)))
  
  stopifnot(selected_transcript %in% rownames(norm_count_matrix) )
  
  transcript_exp<- norm_count_matrix[which(rownames(norm_count_matrix)==selected_transcript),]
  transcript_exp_df<- data.frame(transcript_exp, stringsAsFactors = F)
  head(transcript_exp_df)
  colnames(transcript_exp_df)<- selected_transcript
  
  transcript_exp_df$ID<-rownames(transcript_exp_df)
  
  phenotype$ID<- rownames(phenotype)
  
  pheno_exp<- left_join(phenotype,transcript_exp_df, by="ID" )
  
  
  message(paste0("Generate violin plot for ",selected_transcript," and stratified by ",strata ))
  
  

  
  violin<- ggplot(pheno_exp, aes_string(x = strata, y =selected_transcript ,group = strata, fill = strata)) + # I prefer to store all the aes() in the first ggplot() layer so that the remaining layers can just be about customising the plot
           geom_violin(trim = FALSE,alpha = 0.5, draw_quantiles=c(0.5),position = position_dodge(1)) +
           geom_boxplot(width = 0.1,position = position_dodge(1)) +
          theme_bw()+ theme(legend.position = "none")
  
  
  return(violin)
  
}



#' Title Volcano plot
#'
#' @param deg_res A data frame of results differential of expression genes
#' @param significant_threshold  threshold for significant genes based on empirical-padj bh, example: 0.1
#' @return volcano plot 
#' 
#' @export 
#'
#' @examples
#' library(ggrepel)
#' library(dplyr)
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' data(phenotype)
#' trait<-"Trait.1"
#' covars<-"Age+Sex"
#' 
#' median_norm<- median_normalization(rnaseq_count_matrix)
#' clean_count_matrix <- apply_filters(count_matrix = median_norm, median_min = 1, expression_sum_min = 10, 
#'                                    max_min = 10, range_min = 5, prop_zero_max = 0.5)
#' deg_res<-lm_count_mat_emp_pval(count_matrix=clean_count_matrix, pheno=phenotype, trait=trait, covariates_string=covars, 
#'                              n_permute=100, log_transform = "log_replace_half_min",
#'                              seed = NULL,outcome_type="continuous")
#' volcano_plot(deg_res=deg_res,significant_threshold=0.1 )



volcano_plot<- function(deg_res, significant_threshold){
  
  
  deg_res$expression = ifelse(deg_res$bh_emp_pvals < significant_threshold & abs(deg_res$adjLogFC) >= 0, 
                              ifelse(deg_res$adjLogFC> 0 ,"Up-regulated",'Down-regulated'),"Stable")
  

  sig_line<- max(deg_res$emp_pvals[which(deg_res$bh_emp_pvals < significant_threshold)])
  
  message("Generate volcano plot..")
  
  deg_res$sig<- ifelse(deg_res$bh_emp_pvals<significant_threshold, "annotate", NA) #Will have different colors depending on significance
  deg_res<- deg_res %>% dplyr::arrange(emp_pvals)
  
  volcano <- ggplot(data = deg_res,  aes(x = adjLogFC, y = -log10(emp_pvals),colour=expression)) +
             geom_point(alpha=0.4, size=3.5) +
             scale_color_manual(labels = c("Down-regulated", "Stable","Up-regulated"), values=c("blue", "grey","red"))+
             geom_hline(yintercept=-log10(sig_line),lty=4,col="black",lwd=0.8) +
             labs(x="log2(folds change per 1 unit increase in exposure)",
             y="-log10(empirical p-value)", 
             title="Differential expression")  + theme_bw()+ 
             theme(plot.title = element_text(hjust = 0.5), 
                   legend.position="none", 
                   legend.title = element_blank())+
            geom_label_repel(data=deg_res[!is.na(deg_res$sig),][1:5,], aes(label=as.factor(geneID)), alpha=0.7, size=2.8, force=1.3)
    
  
  return(volcano)
}
