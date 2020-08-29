
#' Title Linear Regressions
#'
#' @param count_matrix : A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#' @param pheno : phenotype data, includes the trait and covariates.  
#' @param trait : A character, the name of the exposure variable. The trait should be a vector in pheno.
#' @param covariates_string : a character string with specifying the covariats, include "as.factor" statements. example: covariate_string = "age + as.factor(sex)"
#' @param log_transform_type  : One of the transformations log_replace_min, log_add_min, log_add_0.5, or NULL (default)
#' @param geneIDs : vector of selection of geneID, NULL if all genes are tested
#' @return Linear regression results as a data frame with columns GeneID, Beta,SE,T_stat (t-statistic),P_value, Z_score (transformed by p-value)
#' @export
#'
#' @examples
#' 

run_linear_regression <-function(count_matrix, pheno, trait, covariates_string, normal_pval = FALSE, geneIDs=NULL, log_transform = NULL){
  
  stopifnot(colnames(count_matrix) == rownames(pheno), is.element(trait, colnames(pheno)))
  
  count_matrix <- as.matrix(count_matrix)
  
  # if geneIDs are (or is) provided, filter count_matrix to the requested genes
  # make sure that requested genes are in the data!
  if(!is.null(geneIDs)){
    if (!all(geneIDs %in% rownames(count_matrix))) {
      genes_not_avail <- geneIDs[which(!geneIDs %in% rownames(count_matrix))]
      message(paste("Some requested geneIDs not in count_matrix, proceeding without using ", 
                        paste(genes_not_avail, collapse = ", ")))
    }
    genes_avail <- geneIDs[which(geneIDs %in% rownames(count_matrix))]
    if (length(genes_avail) == 0){
      stop("None of the requested genes are in the count_matrix, stopping...")
    }
    message(paste0("Filtering count_matrix to genes", genes_avail))
              
    count_matrix<-count_matrix[rownames(count_matrix) %in% geneIDs, , drop = FALSE]
    }
  
  count_matrix <- log_transform_count(count_matrix, transform = log_transform_type)

  # transpose the matrix of counts:
  count_matrix <- t(count_matrix)

  model_string <- paste(covariates_string, "+", trait)
  # 
  XX<-model.matrix(as.formula(paste0("~", model_string)), data=pheno)
  XtXinv <- solve(t(XX) %*% as.matrix(XX))
  XtXinv_se_arg <- sqrt(XtXinv[trait,trait])
  numExplan <-ncol(XX)
  
  XXproj <- XtXinv %*% t(XX)
  betas_mat <- XXproj %*% count_matrix
  betas <- betas_mat[trait,]
  
  resid_Ys <- count_matrix - XX %*% XXproj %*% count_matrix
  sum_squares_resids <- colSums(resid_Ys^2)
  sigmas_square <- sum_squares_resids/(nrow(count_matrix)-numExplan)
  se_betas <- sqrt(sigmas_square)*XtXinv_se_arg
  
  test_stats <- betas/se_betas
  t_pval <- 2*pt(abs(test_stats), lower.tail=FALSE, df = nrow(count_matrix) - numExplan)
  res <- data.frame(outcome = colnames(count_matrix), beta = betas, se_beta = se_betas, test_stat = test_stats, pvalue = t_pval)

  colnames(res)<- c("GeneID", "Beta","SE","T_stat","P_value")
  rownames(res)<-NULL 
  res<- res %>% mutate(FDR_BH= p.adjust(Pvalue, method = "BH"), Z_score= qnorm(1-(P_value/2))*sign(Beta))
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




#' Title Filter Genes
#'
#' @param gene_counts: raw counts (After harmonization based on Phenotype) Ronames gene_counts has to be identical with phnotype ID
#'                       - rowname must be ENSMBL ID :
#'                         example :
#'                         TOR841324 TOR127830 TOR257836 TOR461713 TOR508155
#'      ENSG00000000003        15        20        16         6         9      
#'      ENSG00000000005         0         0         0         0         0
#'      ENSG00000000419       659       832       855       704       564
#'      ENSG00000000457       595       946       721       656       555
#'      
#' @param cv : Covariance variation (SD/mean)
#' @param median_count  : median
#' @param mmr : Max/media ----> inf value will replace to 0
#' @param percent_zero_count  : percentage of 0 allowed in gen count (count_matrix)
#' @param range : max-min 
#' @param pheno  : phenotype after harmonizations
#'
#' @return (filter gene count or matrix counts )

filter_genes<- function(gene_counts, cv_max =NULL, cv_min=NULL, median_val= NULL,mmr_val=NULL, percent_zero_count= NULL,range_val=NULL,Q1_val=NULL, Q3_val=NULL, IQR_val=NULL ,pheno){
  
  
  #normalized the data
  norm_count<- normalized(gene_counts)
  dim(norm_count)
  
  ############ Filter method
  
  #countsx<-count_matrix
  
  zero_count<- rowSums(norm_count == 0)
  
  zero_count<- data.frame(zero_count= zero_count) %>% rownames_to_column(var="gene_id")
  head(zero_count, 10L)
  
  #Median
  median_count<- apply(norm_count, 1, median)
  
  median_count<- data.frame(median_count= median_count) %>%rownames_to_column(var="gene_id")
  
  # Maximum count
  max_count<- apply(norm_count, 1, max)
  max_count<- data.frame(max_count= max_count)%>%rownames_to_column(var="gene_id")
  
  #std deviation count
  sd_count<- apply(norm_count, 1, sd)
  sd_count<- data.frame(sd_count= sd_count)%>%rownames_to_column(var="gene_id")
  
  #Mean count
  mean_count<- apply(norm_count, 1, mean)
  mean_count<- data.frame(mean_count= mean_count)%>%rownames_to_column(var="gene_id")
  
  # Minimum Count
  min_count<- apply(norm_count, 1, min)
  min_count<- data.frame(min_count= min_count)%>%rownames_to_column(var="gene_id")
  
  
  #Q1
  Q1_count<- apply(norm_count, 1, quantile, prob=0.25) 
  Q1_count<- data.frame(Q1= Q1_count)%>%rownames_to_column(var="gene_id")
  
  Q3_count<- apply(norm_count, 1, quantile, prob=0.75) 
  Q3_count<- data.frame(Q3= Q3_count)%>%rownames_to_column(var="gene_id")
  
  IQR_count<-apply(norm_count, 1, IQR) 
  IQR_count<- data.frame(IQR= IQR_count)%>%rownames_to_column(var="gene_id")
  
  
  
  gene_count<- rownames_to_column(as.data.frame(norm_count), var = "gene_id")
  head(gene_count)
  
  #join all for median, mean, max, min, zero count ,and sd
  gene_summary<- join_all(list(gene_count,zero_count,median_count,max_count, sd_count,mean_count,min_count,Q1_count,Q3_count,IQR_count), by="gene_id", type="left")
  dim(gene_summary)
  
  #gene_summary[,464:469]
  
  gene_summary<- gene_summary %>% dplyr::mutate(cv=sd_count/mean_count ) %>%
    dplyr::mutate(mmr=max_count/median_count ) %>%
    dplyr::mutate(range=max_count-min_count)
  
  
  ## This Early filter to remove lowly express (remove all the gene which have max value < 10)
  gene_summary<-gene_summary[which(gene_summary$max_count > 10),]
  dim(gene_summary)
  
  if(!is.null(cv_max)){
    #remove median below 5
    counts_cv_max<-gene_summary[which(gene_summary$cv <= cv_max),]
    dim(counts_cv_max)
    gene_summary<- counts_cv_max
    
  }
  
  if(!is.null(cv_min)){
    #remove median below 5
    counts_cv_min<-gene_summary[which(gene_summary$cv >= cv_min),]
    dim(counts_cv_min)
    gene_summary<- counts_cv_min
    
  }
  
  if(!is.null(percent_zero_count)){
    
    counts_nonzero<-gene_summary[which(gene_summary$zero_count<= percent_zero_count*nrow(pheno)),]
    gene_summary<- counts_nonzero
    #
  }
  #Filter Median
  if(!is.null(median_val)){
    counts_med<-gene_summary[which(gene_summary$median_count >= median_val),]
    dim(counts_med)
    gene_summary<- counts_med
    
  }
  ## MMR
  if(!is.null(mmr_val)){
    gene_summary$mmr[is.infinite(gene_summary$mmr)] <- 0
    counts_mmr<- gene_summary[which(gene_summary$mmr < mmr_val),]
    gene_summary<- counts_mmr
    
  }
  ### range
  if(!is.null(range_val)){
    counts_range<- gene_summary[which(gene_summary$range < range_val),]
    gene_summary<- counts_range
    
  }
  
  if(!is.null(Q1_val)){
    counts_Q1<- gene_summary[which(gene_summary$Q1 > Q1_val ),]
    gene_summary<- counts_Q1
    
  }
  
  if(!is.null(Q3_val)){
    counts_Q3<- gene_summary[which(gene_summary$Q3 < Q3_val ),]
    gene_summary<- counts_Q3
    
  }
  
  if(!is.null(IQR_val)){
    counts_IQR<- gene_summary[which(gene_summary$IQR < IQR_val ),]
    gene_summary<- counts_IQR
    
  }
  
  print(dim(gene_summary))
  
  
  counts<- norm_count[gene_summary$gene_id,]
  return(counts)
}

