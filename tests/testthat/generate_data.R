.testData<- function(n_gene=1000, n_samples=20){
  set.seed(0)
  m0 <- matrix(0, n_gene, n_samples)
  
  #count_mat <- matrix(as.integer(rexp(20000, rate=.01)), ncol=20)
  
  
  #dim(b)
  
  count_mat<-apply(m0, c(1,2), function(x) sample(c(0,1:n_gene),1))

  rownames(count_mat)<- paste0("gene_",1:n_gene)
  colnames(count_mat)<- paste0("sample_",1:n_samples)
  
  
  pheno<- data.frame(Age=runif(n_samples,40,100),
                     Sex=sample(c("F","M"), n_samples, replace=TRUE, prob=c(0.7, 0.3)),
                     Trait.1=rnorm(n_samples,10,5), 
                     Trait.2= rnorm(n_samples,4,2), 
                     Trait.3=sample(c(1,0), n_samples, replace=TRUE, prob=c(0.6, 0.4)))
  
  rownames(pheno)<- colnames(count_mat)
  
  return(list(gene_count=count_mat, phenotype=pheno))
}


