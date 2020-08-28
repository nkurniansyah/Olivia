library(BatchJobs)
library(parallel)
library(GWASTools)
library(dplyr)
library(tidyverse)
library(plyr)
library(foreach)
library(doParallel)
library(progress)
library(qvalue)
library(EnsDb.Hsapiens.v86)
library(argparser)
library(data.table)
library(gtools)
library(EnsDb.Mmusculus.v79)




source("/Volumes/linkage/MESA2/Projects/2019-RNA-seq/linear-regression_rna-seq/gene_filters.R")
source("/Volumes/linkage/MESA2/Projects/2019-RNA-seq/linear-regression_rna-seq/Utils.R")
source("/Volumes/linkage/MESA2/Projects/2019-RNA-seq/linear-regression_rna-seq/method.R")

argp <- arg_parser("Run LR with Residual Normalizations")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--permute", help="seed (number of permutations)", type="integer")

argv <- parse_args(argp)
config <- readConfig(argv$config)
permute <- argv$permute

required <- c("trait",
              "pheno_file",
              "residual_perm_path",
              "output",
              "gene_counts_file",
              "covariates_string")
optional <- c("median_val"= NA, # takes precedence
              "cv_max"=NA,
              "cv_min"=NA,
              "mmr_val"= NA,
              "percent_zero_count"=NA,
              "range_val"=NA,
              "family"="gaussian",
              "list_geneID"=NA,
              "Q1_val"= NA,
              "Q3_val"=NA,
              "IQR_val"=NA,
              "parallel"="TRUE",
              "organism"="Mmusculus")


config <- setConfigDefaults(config, required, optional)
print(config)

parallel<- config["parallel"]

trait <- as.character(config["trait"])

start.time <- Sys.time()


#pheno<- read.csv("/Volumes/linkage/MESA2/Projects/2019-RNA-seq/RNAseq/Data/2020-08-28_pheno.csv")
pheno<- read.csv(config["pheno_file"])
head(pheno)
#pheno<- pheno %>% dplyr::rename(ID=TOR_ID)


#mat<- fread("/Volumes/linkage/MESA2/Projects/2019-RNA-seq/RNAseq/Data/GSE156437_RNA_seq.txt", data.table = F)
## in RData



mat<- fread(config["gene_counts_file"], data.table = F)

#mat[1:10,1:10]

row.names(mat)<- mat[,1]

mat<-mat[,-1]

mat<-as.matrix(mat)

storage.mode(mat)<-"numeric"
##### Fix Your matrix format in here
#rownames(mat)<- mat[,1]
#mat<- mat[,-1]

#dim(mat)

organism<-config["organism"]
median_val<- as.numeric(config["median_val"])
cv_max<- as.numeric(config["cv_max"])
cv_min<- as.numeric(config["cv_min"])
mmr_val<- as.numeric(config["mmr_val"])
percent_zero_count<-as.numeric(config["percent_zero_count"])
range_val<-as.numeric(config["range_val"])
Q1_val<- as.numeric(config["Q1_val"])
Q3_val<- as.numeric(config["Q3_val"])
IQR_val<- as.numeric(config["IQR_val"])

cv_max<- NA
cv_min<- NA
median_val<- NA
mmr_val<-NA
percent_zero_count<-0.5
range_val<-NA
Q1_val<- NA
Q3_val<-NA
IQR_val<-NA

if(is.na(cv_max)){
  cv_max<-NULL
}else {
  cv_max<- as.numeric(cv_max)
}

if(is.na(cv_min)){
  cv_min<-NULL
}else {
  cv_min<- as.numeric(cv_min)
}
if(is.na(median_val)){
  median_val<-NULL
}else {
  median_val<- as.numeric(median_val)
}

if(is.na(mmr_val)){
  mmr_val<-NULL
}else {
  mmr_val<- as.numeric(mmr_val)
}

if(is.na(percent_zero_count)){
  percent_zero_count<-NULL
}else {
  percent_zero_count<- as.numeric(percent_zero_count)
}

if(is.na(range_val)){
  range_val<-NULL
}else {
  range_val<- as.numeric(range_val)
}

if(is.na(Q1_val)){
  Q1_val<-NULL
}else {
  Q1_val<- as.numeric(Q1_val)
}

if(is.na(Q3_val)){
  Q3_val<-NULL
}else {
  Q3_val<- as.numeric(Q3_val)
}

if(is.na(IQR_val)){
  IQR_val<-NULL
}else {
  IQR_val<- as.numeric(IQR_val)
}

genes_filter_mat<- filter_genes(gene_counts = mat , cv_max = cv_max, cv_min=cv_min, median_val = median_val,mmr_val = mmr_val, percent_zero_count = percent_zero_count, range_val = range_val,Q1_val=Q1_val, Q3_val=Q3_val, IQR_val=IQR_val, pheno = pheno )
dim(genes_filter_mat)

#head(pheno)

#no space in covariantes string: Age,BMI,Sex

covariates_string<- as.character(config["covariates_string"])
print(covariates_string)

covars<- gsub(",","+",covariates_string)



list_geneID<- config["list_geneID"] # If you want to run only specific transcripts 
family <- as.character(config["family"]) # value of the family trait  ( Gaussian or binomial)

if(!is.na(list_geneID)){
  list_geneID<- read.csv(list_geneID, header = F)
  list_geneID<- as.vector(list_geneID$V1)
  results<- run_linear_regression(count_matrix = genes_filter_mat, pheno =pheno,trait = trait ,covariates_string = covariates_string,log_transform = log_replace_min, list_geneID = list_geneID)
  head(results) 
  gene_symbol<- select(EnsDb.Hsapiens.v86, keys =as.character(results$GeneID) , keytype = "GENEID",
                       columns = c("GENEID", "GENENAME"))
  
  head(gene_symbol)
  
  colnames(gene_symbol)<- c("GeneID","GeneName")
  
  
  #### Repalce unknown gene symbel with ENSG ID
  ind<-which(is.na(gene_symbol$GeneName))
  if(any(is.na(gene_symbol$GeneName))){
    temp_gene_name<- gene_symbol$GeneID
    gene_symbol$GeneName[ind] <- temp_gene_name[ind]
  }else{
    gene_symbol<-gene_symbol
  }
  
  annot_genes<- left_join(results,gene_symbol ,by="GeneID" )
  head(annot_genes)
  
  annot_genes<- annot_genes %>% dplyr::select(GeneName, Beta,SE,Stat,z.score,Pvalue)
  dim(annot_genes)
  head(annot_genes)
  
  write.csv(annot_genes, file= config["output"], row.names = F)
  
  
} else {
  list_geneID <- NULL
  #setup parallel backend to use many processors
  
  results<- run_linear_regression(count_matrix = genes_filter_mat, pheno =pheno,trait = trait ,covariates_string = covariates_string,log_transform = log_replace_min, list_geneID = list_geneID)
  head(results) 
  
  
  
  results[which(results$FDR_BH < 0.05),]
  
  if (parallel== "TRUE"){
    cores<-detectCores()
    tocor<- cores-1
    print(paste0("Using ",tocor," core"))
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    
    registerDoParallel(cl)
    print(cl)
    
    
    
    
    # progress bar ------------------------------------------------------------
    
    iterations <- 10
    pb <- txtProgressBar(0, iterations, style = 2)
    
    list<-foreach(permute=1:iterations) %do% {
      start<- ((permute-1) * 10) +1
      end<- (permute) * 10
      setTxtProgressBar(pb, permute)
      for(i in start:end){
        
        if(family=="gaussian"){
          perm_res<- permute_resids_trait(pheno=pheno, outcome = trait ,covariates_string =covars, seed = i)
          
        }else if (family=="binomial"){
          perm_res<- permute_resids_trait_bino(pheno=pheno, outcome = trait ,covariates_string =covars, seed = i, family = "binomial")
        }
        phenotype<- cbind(pheno, perm_res)
        
        residual_permutations<-run_linear_regression(count_matrix =genes_filter_mat,covariates_string = covariates_string, trait = "perm_res" ,pheno = phenotype,log_transform = log_replace_min)
        head(residual_permutations)
        residual_perm_dir<- "/Volumes/linkage/MESA2/Projects/2019-RNA-seq/RNAseq/Residual_Permutation"
        #print(paste0(output_path,"/Residual_perm_",trait,"_",i,".RData"))
        #residual_perm_dir<- config["residual_perm_path"]
        save(residual_permutations, file = paste0(residual_perm_dir,"/Residual_perm_",trait,"_",i,".RData"), version = 2)
      }
      
    }
    
    stopCluster(cl)
    
  }else{
    
    print("no paralel")
    
    for(i in 1:100){
      
      if(family=="gaussian"){
        perm_res<- permute_resids_trait(pheno=pheno, outcome = trait ,covariates_string =covars, seed = i)
        
      }else if (family=="binomial"){
        perm_res<- permute_resids_trait_bino(pheno=pheno, outcome = trait ,covariates_string =covars, seed = i, family = "binomial")
      }
      phenotype<- cbind(pheno, perm_res)
      
      residual_permutations<-run_linear_regression(count_matrix =genes_filter_mat,covariates_string = covariates_string, trait = "perm_res" ,pheno = phenotype,log_transform = log_replace_min)
      #output_path<- "/data/linkage/MESA2/Projects/2019-RNA-seq/Analysis/PBMC/Example/Residual_perm"
      #print(paste0(output_path,"/Residual_perm_",trait,"_",i,".RData"))
      residual_perm_dir<- config["residual_perm_path"]
      save(residual_permutations, file = paste0(residual_perm_dir,"/Residual_perm_",trait,"_",i,".RData"))
    }
    
  }
  
  
  
  ###### compute Quantile empirical Pval
  
  null_pval_files<- dir(residual_perm_dir, full.names = T, recursive = F, pattern = ".RData")
  null_pval_files<- mixedsort(null_pval_files)
  
  stopifnot(length(null_pval_files) <= 100)
  
  null_pvals <-  vector("list", nrow(results)*100)
  
  for(i in 1:100){
    results_id<- results %>% dplyr::select(GeneID)
    #print(null_pval_files[i])
    perm_pvals_df<- getobj(null_pval_files[i])
    perm_pvals_df<- left_join(results_id,perm_pvals_df, by="GeneID")
    perm_pvals_df<-perm_pvals_df$Pvalue
    null_pvals[[i]]<-perm_pvals_df
  }
  null_pval_res<- unlist(null_pvals)
  length(null_pval_res)
  
  # get emperical pval
  emppvals<- quantile_emPval(p_values = results$Pvalue, null_p_values = null_pval_res)
  
  
  emp_pvals_df<- cbind(results, emPval=emppvals)
  head(emp_pvals_df)
  
  #add emp FDR BH, FWER home , qval and lfdr
  emp_pvals_df<- emp_pvals_df %>% mutate(emp_FDR_BH=p.adjust(emPval, method = "BH")) %>%
    mutate(emp_fwer_holm=p.adjust(emPval, method = "holm"))%>%
    mutate(emp_lfdr=qvalue(emPval)$lfdr)
  
  dim(emp_pvals_df)
  head(emp_pvals_df)
  
  results$FDR_BH
  emp_pvals_df[which(emp_pvals_df$emp_FDR_BH < 0.05),]
  
  #edb <- EnsDb.Hsapiens.v86
  
  
  #unknown_genes<- deg_avgsat %>% dplyr::filter(str_detect(gene_names,"ENS"))
  if(organism=="Mmusculus"){
    orgDB<-EnsDb.Mmusculus.v79
  }else if(organism=="Hsapiens"){
    orgDB<-EnsDb.Hsapiens.v86
  }
  
  gene_symbol<- select(orgDB, keys =as.character(emp_pvals_df$GeneID) , keytype = "GENEID",
                       columns = c("GENEID", "GENENAME"))
  
  head(gene_symbol)
  
  colnames(gene_symbol)<- c("GeneID","GeneName")
  
  
  #### Repalce unknown gene symbel with ENSG ID
  ind<-which(is.na(gene_symbol$GeneName))
  if(any(is.na(gene_symbol$GeneName))){
    temp_gene_name<- gene_symbol$GeneID
    gene_symbol$GeneName[ind] <- temp_gene_name[ind]
  }else{
    gene_symbol<-gene_symbol
  }
  
  annot_genes<- left_join(emp_pvals_df,gene_symbol ,by="GeneID" )
  head(annot_genes)
  
  annot_genes<- annot_genes %>% dplyr::select(GeneName, Beta,SE,Stat,z.score,Pvalue,FDR_BH,emPval ,emp_FDR_BH)
  dim(annot_genes)
  head(annot_genes)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  #remove Resid permute files
  
  file.remove(null_pval_files)
  
  #out[[i]]<- add_val
  write.csv(annot_genes, file= config["output"], row.names = F)
  #hist(emppvals)
  
}
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")

