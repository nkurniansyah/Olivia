library(BatchJobs)
library(parallel)
library(ggplot2)
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

source("/data/linkage/MESA2/Projects/2019-RNA-seq/linear-regression_rna-seq/gene_filters.R")
source("/data/linkage/MESA2/Projects/2019-RNA-seq/linear-regression_rna-seq/Utils.R")
source("/data/linkage/MESA2/Projects/2019-RNA-seq/linear-regression_rna-seq/method.R")


argp <- arg_parser("Run LR with Residual Normalizations")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--permute", help="seed (number of permutations)", type="integer")

argv <- parse_args(argp)
config <- readConfig(argv$config)
permute <- argv$permute

required <- c("pheno_file",
              "residual_perm_path",
              "output",
              "gene_counts_file",
              "covariates_string",
              "exposures")
optional <- c("median_val"= NA, # takes precedence
              "cv_max"=NA,
              "cv_min"=NA,
              "mmr_val"= NA,
              "percent_zero_count"=NA,
              "range_val"=NA,
              "family"="gaussian",
              "list_geneID"=NA,
              "parallel"="TRUE",
              "Q1_val"= NA,
              "Q3_val"=NA,
              "IQR_val"=NA,
              "organism"="Mmusculus")


config <- setConfigDefaults(config, required, optional)
print(config)
parallel<- config["parallel"]


start.time <- Sys.time()
#pheno<- read.csv("/data/linkage/MESA2/Projects/2019-RNA-seq/Analysis/PBMC/Data/2020-03-09_example_pheno.csv")
pheno<- read.csv(config["pheno_file"])
head(pheno)
#pheno<- pheno %>% dplyr::rename(ID=TOR_ID)
mat<- fread(config["gene_counts_file"], data.table = F)

#mat[1:10,1:10]

row.names(mat)<- mat[,1]

mat<-mat[,-1]

mat<-as.matrix(mat)

storage.mode(mat)<-"numeric"

organism<-config["organism"]

if(organism=="Mmusculus"){
  orgDB<-EnsDb.Mmusculus.v79
}else if(organism=="Hsapiens"){
  orgDB<-EnsDb.Hsapiens.v86
}



median_val<- config["median_val"]
cv_max<- config["cv_max"]
cv_min<- config["cv_min"]
mmr_val<- config["mmr_val"]
percent_zero_count=config["percent_zero_count"]
range_val= config["range_val"]
Q1_val<- config["Q1_val"]
Q3_val<- config["Q3_val"]
IQR_val<- config["IQR_val"]


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
genes_filter_mat[1:10,1:10]

head(pheno)

#covariates_string NO SPACE
covariates_string<- as.character(config["covariates_string"])
print(covariates_string)

covars<- gsub(",","+",covariates_string)

exposures<- config["exposures"]
print(exposures)



list_geneID<- config["list_geneID"]

family <- as.character(config["family"])
if(!is.na(list_geneID)){
  list_geneID<- read.csv(list_geneID, header = F)
  list_geneID<- as.vector(list_geneID$V1)
  
  exposures<- strsplit(exposures,",")[[1]]
  
  results<- run_multivar_Linear_reggressions(count_matrix = genes_filter_mat, pheno =pheno,exposures = exposures ,covariates_string = covariates_string, log_transform = log_replace_min, list_geneID = list_geneID)
  head(results) 
  gene_symbol<- select(orgDB, keys =as.character(results$GeneID) , keytype = "GENEID",
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
  
  annot_genes<- annot_genes %>% dplyr::select(GeneName, contains("Beta"),Joint_stats,Joint.Pvalue)
  dim(annot_genes)
  head(annot_genes)
  
  write.csv(annot_genes, file= config["output"], row.names = F)
  
  
} else {
  list_geneID <- NULL
  #setup parallel backend to use many processors
  exposures<- strsplit(exposures,",")[[1]]
  
  results<- run_multivar_Linear_reggressions(count_matrix = genes_filter_mat, pheno =pheno,exposures = exposures ,covariates_string = covariates_string, log_transform = log_replace_min, list_geneID = list_geneID)
  head(results) 
  
  if (parallel== "TRUE"){
    
    
    cores=detectCores()
    tocor<- cores-1
    print(paste0("Using ",tocor," core"))
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
    
    
    # progress bar ------------------------------------------------------------
    
    iterations <- 10
    pb <- txtProgressBar(0, iterations, style = 2)
    
    list<-foreach(permute=1:iterations) %do% {
      
      start<- ((permute-1) * 10) +1
      end<- (permute) * 10
      setTxtProgressBar(pb, permute)
      for(i in start:end){
        exposures<- gsub(","," ",exposures)
        exposures<- unlist(strsplit(exposures,split = " "))
        
        exp<- list()
        for (exposure in exposures){
          
          if(family=="gaussian"){
            
            perm_res<- permute_resids_trait(pheno=pheno, outcome = exposure ,covariates_string =covars, seed = i)
            
          }else if (family=="binomial"){
            perm_res<- permute_resids_trait_bino(pheno=pheno, outcome = exposure ,covariates_string =covars, seed = i, family = "binomial")
          }
          perm_res<- data.frame(perm_res)
          colnames(perm_res)<- paste0("perm_",exposure)
          exp[[exposure]]<-perm_res
        }
        expo<- do.call(cbind,exp)
        phenotype<- cbind(pheno, expo)
        
        resid_exposure<- paste0("perm_",exposures)
        residual_permutations<- run_multivar_Linear_reggressions(count_matrix = genes_filter_mat, pheno =phenotype,exposures = resid_exposure ,covariates_string = covariates_string, log_transform = log_replace_min, list_geneID = list_geneID)
        
        print(head(residual_permutations))
        #output_path<- "/data/linkage/MESA2/Projects/2019-RNA-seq/Analysis/PBMC/Example/Residual_perm"
        #print(paste0(output_path,"/Residual_perm_",trait,"_",i,".RData"))
        residual_perm_dir<- config["residual_perm_path"]
        save(residual_permutations, file = paste0(residual_perm_dir,"/Residual_perm_multivar_",i,".RData"))
      }
      
    }
    
    stopCluster(cl)
    
  }else{
    print("no paralel")
    
    for(i in 1:100){
      exposures<- gsub(","," ",exposures)
      exposures<- unlist(strsplit(exposures,split = " "))
      
      exp<- list()
      for (exposure in exposures){
        
        if(family=="gaussian"){
          
          perm_res<- permute_resids_trait(pheno=pheno, outcome = exposure ,covariates_string =covars, seed = i)
          
        }else if (family=="binomial"){
          perm_res<- permute_resids_trait_bino(pheno=pheno, outcome = exposure ,covariates_string =covars, seed = i, family = "binomial")
        }
        perm_res<- data.frame(perm_res)
        colnames(perm_res)<- paste0("perm_",exposure)
        exp[[exposure]]<-perm_res
      }
      expo<- do.call(cbind,exp)
      phenotype<- cbind(pheno, expo)
      
      resid_exposure<- paste0("perm_",exposures)
      residual_permutations<-run_multivar_Linear_reggressions(count_matrix =genes_filter_mat,covariates_string = covariates_string, exposures = resid_exposure ,pheno = phenotype,log_transform = log_replace_min, list_geneID =list_geneID )
      head(residual_permutations)
      #output_path<- "/data/linkage/MESA2/Projects/2019-RNA-seq/Analysis/PBMC/Example/Residual_perm"
      #print(paste0(output_path,"/Residual_perm_",trait,"_",i,".RData"))
      residual_perm_dir<- config["residual_perm_path"]
      save(residual_permutations, file = paste0(residual_perm_dir,"/Residual_perm_multivar_",i,".RData"))
    }
  }
  
  
  ###### compute Quantile empirical Pval
  
  null_pval_files<- dir(residual_perm_dir, full.names = T, recursive = F, pattern = ".RData")
  null_pval_files<- mixedsort(null_pval_files)
  
  stopifnot(length(null_pval_files) <= 100)
  
  null_pvals <-  vector("list", nrow(results)*100)
  
  for(i in 1:100){
    results_id<- results %>% dplyr::select(GeneID)
    print(null_pval_files[i])
    perm_pvals_df<- getobj(null_pval_files[i])
    perm_pvals_df<- left_join(results_id,perm_pvals_df, by="GeneID")
    perm_pvals_df<-perm_pvals_df$Joint.Pvalue
    null_pvals[[i]]<-perm_pvals_df
  }
  null_pval_res<- unlist(null_pvals)
  length(null_pval_res)
  
  # get emperical pval
  emppvals<- quantile_emPval(p_values = results$Joint.Pvalue, null_p_values = null_pval_res)
  
  
  emp_pvals_df<- cbind(results, Joint.emPval=emppvals)
  head(emp_pvals_df)
  
  #add emp FDR BH, FWER home , qval and lfdr
  emp_pvals_df<- emp_pvals_df %>% mutate(Joint.emp_FDR_BH=p.adjust(Joint.emPval, method = "BH")) %>%
    mutate(Joint.emp_fwer_holm=p.adjust(Joint.emPval, method = "holm"))%>%
    mutate(Joint.emp_lfdr=qvalue(Joint.emPval)$lfdr)
  
  dim(emp_pvals_df)
  head(emp_pvals_df)
  
  #edb <- EnsDb.Hsapiens.v86
  
  
  #unknown_genes<- deg_avgsat %>% dplyr::filter(str_detect(gene_names,"ENS"))
  
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
  
  annot_genes<- annot_genes %>% dplyr::select(GeneName, contains("Beta"),Joint_stats,Joint.Pvalue, Joint.FDR_BH,Joint.emPval,Joint.emp_FDR_BH)
  dim(annot_genes)
  head(annot_genes)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  
  #out[[i]]<- add_val
  write.csv(annot_genes, file= config["output"], row.names = F)
  #hist(emppvals)
  
}

ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")

