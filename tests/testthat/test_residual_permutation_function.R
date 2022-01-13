


test_that("Residual permutation function trait continous", {
  set.seed(123)
  data<-.testData()
  pheno<- data$phenotype
  trait<-"Trait.1"
  covariates<-"Age+Sex" # change to covariates_string which is like a formula
  
  
  perm_resid<-permute_resids_trait(pheno = pheno,trait = trait,
                                   seed = NULL,covariates_string =covariates,
                                   outcome_type  = "continuous")
  
  
  expect_true(all(names(perm_resid) %in% rownames(pheno)), TRUE)
  
  pheno$binary_out <- rbinom(nrow(pheno), 1, pheno$Age/max(pheno$Age))
  sampled_binary_out <- permute_resids_trait(pheno = pheno,trait = "binary_out",
                                   seed = NULL,covariates_string =covariates,
                                   outcome_type  = "binary")
  
  expect_true(all(sampled_binary_out %in% c(0,1)), TRUE)
})




test_that("Residual permutation function trait binary", {
  set.seed(123)
  data<-.testData()
  pheno<- data$phenotype
  trait<-"Trait.3"
  covariates<-"Age+Sex" # change to covariates_string which is like a formula
  
  perm_resid<-permute_resids_trait(pheno = pheno,trait = trait,
                                   seed = NULL,covariates_string =covariates,
                                   outcome_type = "binary")
  
  expect_true(all(names(perm_resid) %in% rownames(pheno)), TRUE)
  
  
})



test_that("Residual permutation function trait power", {
  set.seed(123)
  data<-.testData()
  pheno<- data$phenotype
  trait<-"Trait.1"
  covariates<-"Age+Sex"
  
  gene_count<-data$gene_count
  single_gene_exp<-gene_count[sample(1:nrow(gene_count),1),]
  
  
  perm_resid_power<-permute_resids_trait_cor(pheno = pheno,trait = trait,seed = 123,covariates_string =covariates,
                                             outcome_type = "continous",gene_exp =single_gene_exp, required_cor = 0.3 )
  
  expect_true(all(names(perm_resid_power) %in% rownames(pheno)), TRUE)
  

  
})


