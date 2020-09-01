


#' Title Filter transcripts by applying various filters
#'
#' @param count_matrix: raw counts (After harmonization based on Phenotype) Rownames count_matrix has to be identical with phnotype ID
#'                       - rownames must be ENSMBL ID :
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

apply_filters <- function(count_matrix, cv_max =NULL, cv_min=NULL, median_val= NULL,mmr_val=NULL, percent_zero_count= NULL,range_val=NULL,Q1_val=NULL, Q3_val=NULL, IQR_val=NULL ,pheno){
  
  
  #normalize the data
  # TS: which function is this? also, a function should be called "normalize", not "normalized". 
  norm_count<- normalized(count_matrix) 
  
  ## all parameter should be computed and only then filters should be applied. Because applying
  ## a single filter will change the values computed for the next filter.
  
  
  ############ Compute characteristics of the normalized count matrix. 
  ## TS: why of the normalized? is that what we decided? 
  
  
  zero_count<- rowSums(norm_count == 0)
  
  zero_count<- data.frame(zero_count= zero_count) %>% rownames_to_column(var="gene_id")
  
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

