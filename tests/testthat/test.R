.testData<- function(row=1000, col=20){
  set.seed(0)
  m0 <- matrix(0, row, col)

  count_mat<-apply(m0, c(1,2), function(x) sample(c(0,1:row),1))

  rownames(count_mat)<- paste0("gene.",1:row)
  colnames(count_mat)<- paste0("sample.",1:col)

  count_mat[1:10,1:10]


  pheno<- data.frame(Age=runif(col,40,100),Sex=sample(c("F","M"), col, replace=TRUE, prob=c(0.7, 0.3)),
                     Trait.1=rnorm(col,10,5), Trait.2= rnorm(col,4,2), Trait.3=sample(c("1","0"), col, replace=TRUE, prob=c(0.6, 0.4)))

  rownames(pheno)<- colnames(count_mat)

  return(list(gene_count=count_mat, phenotype=pheno))
}




test_that("filter by gene", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count
  set.seed(5); genes<- sample(rownames(gene_count), 5)
  expect_equal(length(genes), 5)
  expect_true(all(genes  %in% rownames(gene_count)))

  gene_selected<- filter_by_genes(count_matrix = gene_count, gene_IDs =genes )

  expect_equal(nrow(gene_selected),5)

  expect_equal(rownames(gene_selected),c("gene.207","gene.697","gene.715","gene.834","gene.875"))

  expect_message(filter_by_genes(count_matrix = gene_count, gene_IDs =genes ), "Filtering count_matrix to genes : gene.834 gene.875 gene.697 gene.207 gene.715")

})



test_that("Median Normalization", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count

  expect_equal(nrow(gene_count),1000)
  expect_equal(ncol(gene_count),20)
  med_norm<- median_normalization(gene_count)
  expect_true(is.numeric(med_norm), TRUE)
  expect_equal(rowSums(med_norm)[1],c(gene.1 = 8669.21637411038))
  expect_equal(colSums(med_norm)[1],c(sample.1 = 498501.5))

})

#(pheno, trait, covariates_string, seed = NULL, family="gaussian"


test_that("Residual permutation function trait", {
  set.seed(123)
  data<-.testData()
  pheno<- data$phenotype
  trait<-"Trait.1"
  covariates<-"Age,Sex"

  perm_resid<-permute_resids_trait(pheno = pheno,trait = trait,seed = NULL,covariates_string =covariates,
                                   family = "gaussian")

  expect_equal(length(perm_resid), nrow(pheno))

  expect_error(permute_resids_trait(pheno = pheno,trait = trait,seed = NULL,covariates_string =covariates,
                                    family = "binomial"))

  binary_trait<-"Trait.3"

  expect_true(all(na.omit(pheno[,binary_trait]) %in% 0:1))

  expect_error(permute_resids_trait(pheno = pheno,trait = trait,seed = NULL,covariates_string =covariates,
              family = "poisson"),"Requested family type is poisson allowed values are gaussian and binomial " )

  expect_error(permute_resids_trait(pheno = pheno,trait = binary_trait,seed = NULL,covariates_string =("Age,sex"),
                                    family = "binomial"),"covariates not found in the phenotype")


})




test_that("Residual permutation function trait power", {
  set.seed(123)
  data<-.testData()
  pheno<- data$phenotype
  trait<-"Trait.1"
  covariates<-"Age,Sex"

  gene_count<-data$gene_count
  single_gene_exp<-gene_count[sample(1:nrow(gene_count),1),]


  perm_resid_power<-permute_resids_trait_cor(pheno = pheno,trait = trait,seed = 123,covariates_string =covariates,
                                   family = "gaussian",gene_exp =single_gene_exp, required_cor = 0.3 )

  expect_equal(length(perm_resid_power), nrow(pheno))

  expect_error(permute_resids_trait_cor(pheno = pheno,trait = trait,seed = NULL,covariates_string =covariates,
                                    family = "binomial",gene_exp =single_gene_exp, required_cor = 0.3))

  binary_trait<-"Trait.3"

  expect_true(all(na.omit(pheno[,binary_trait]) %in% 0:1))


  expect_error(permute_resids_trait_cor(pheno = pheno,trait = trait,seed = NULL,covariates_string =covariates,
                                      family = "poisson",gene_exp =single_gene_exp, required_cor = 0.3),"Requested family type is poisson allowed values are gaussian and binomial " )

  expect_error(permute_resids_trait_cor(pheno = pheno,trait = binary_trait,seed = NULL,covariates_string =("Age,sex"),
                                    family = "binomial",gene_exp =single_gene_exp, required_cor = 0.3),"covariates not found in the phenotype")


})




test_that("Log replace half min", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count

  gene_count<-gene_count[which(rowSums(gene_count)>0),]
  gene_count_norm<- median_normalization(gene_count)

  log_gene_count<- log_replace_half_min(gene_count)

  expect_equal(dim(log_gene_count),c(1000L, 20L))

  expect_equal(rowSums(log_gene_count)[2],c(gene.2 = 175.107019470175))


  })



test_that("Log replace half min", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count

  gene_count<-gene_count[which(rowSums(gene_count)>0),]
  gene_count_norm<- median_normalization(gene_count)

  log_gene_count<- log_replace_half_min(gene_count)

  expect_equal(dim(log_gene_count),c(1000L, 20L))

  expect_equal(rowSums(log_gene_count)[2],c(gene.2 = 175.107019470175))


})




test_that("Log add half", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count

  gene_count<-gene_count[which(rowSums(gene_count)>0),]
  gene_count_norm<- median_normalization(gene_count)

  log_gene_count<- log_add_0.5(gene_count)

  expect_equal(dim(log_gene_count),c(1000L, 20L))

  expect_equal(rowSums(log_gene_count)[100],c(gene.100 = 166.061984466392))


})


test_that("log Transform Count", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count

  gene_count<-gene_count[which(rowSums(gene_count)>0),]


  gene_count_norm<- log_transform_count(count_matrix =gene_count, transform = "log_replace_half_min")
  expect_equal(dim(gene_count_norm),c(1000L, 20L))

  expect_message(log_transform_count(count_matrix =gene_count, transform = NULL),"No transformation of gene counts requested")

  expect_error(log_transform_count(count_matrix =gene_count, transform = "median"),"Requested transformation not allowed.
         Allowed transformation names are log_replace_half_min, log_add_min, log_add_0.5")


  expect_equal(rowSums(gene_count_norm)[255],c(gene.255 = 176.420952484513))


})


#apply_filters <- function(count_matrix, median_min = 1, expression_sum_min = 10,
#                          max_min = 5, range_min = 5, prop_zero_max = 0.8,
##                          cv_min = NULL, cv_max = NULL,
#                          max_to_median_max = NULL)



test_that("Apply filters", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count

  gene_count<-gene_count[which(rowSums(gene_count)>0),]


  gene_count_norm<- log_transform_count(count_matrix =gene_count, transform = "log_replace_half_min")

  filter_gene<- apply_filters(count_matrix=gene_count_norm, median_min = 1, expression_sum_min = 10,
                            max_min = 5, range_min = 5, prop_zero_max = 0.8,
                            cv_min = NULL, cv_max = NULL,
                            max_to_median_max = NULL)
  expect_equal(dim(filter_gene),c(428L, 20L))

  expect_equal(rowSums(filter_gene)[155],c(gene.361 = 173.40759606542))


})


test_that("Quantile empirical p-values", {
  set.seed(123)
  stat<- runif(1000,0,1)
  null_stat<- runif(100000,0,1)
  emp_pval<-compute_quantile_empirical_pvalues(statistics=stat,null_statistics=null_stat)
  expect_equal(length(stat),length(emp_pval))
  expect_false( stat[2] >emp_pval[2])
  expect_equal(emp_pval[10],0.45884)
})



test_that("Storey empirical p-values", {
  set.seed(123)
  library(qvalue)
  stat<- rnorm(1000,0,1)
  null_stat<- rnorm(100000,0,1)
  sto_emp_pval<-compute_storey_empirical_pvalues(statistics=stat,null_statistics=null_stat)
  expect_equal(length(stat),length(sto_emp_pval))
  expect_false( stat[2] >sto_emp_pval[2])
  expect_equal(sto_emp_pval[10],0.67173)
})



test_that(" Permutation p-values", {
  set.seed(123)
  library(qvalue)
  pval<- 1e-05
  null_pval<-runif(100000,0,1)
  perm_pval<-permutataion_pvalues(pvalue=pval,null_pval=null_pval)
  expect_equal(length(pval),length(perm_pval))
  expect_false( pval < perm_pval)
})



test_that("Computed emprical p-values", {
  set.seed(123)
  library(qvalue)
  stat<- runif(1000,0,1)
  null_stat<- runif(100000,0,1)
  empval<-compute_empirical_pvalues(statistics=stat,null_statistics=null_stat,empirical_type = "quantile",
                            stat_type = "p_value",
                            t_df = NULL)

  expect_message(compute_empirical_pvalues(statistics=stat,null_statistics=null_stat,empirical_type = "quantile",
                                           stat_type = "p_value",
                                           t_df = NULL), "Run quantile empirical p-values using p_value")

  expect_equal(length(stat),length(empval))
  stat_storey<- rnorm(1000,0,3)
  null_stat_storey<- rnorm(100000,0,3)

  expect_message(compute_empirical_pvalues(statistics=stat_storey,null_statistics=null_stat_storey,empirical_type = "storey",
                                           stat_type = "z_score",
                                           t_df = NULL), "Run storey empirical p-values using z_score")

})


test_that("fast Linear regression for single exposure ", {
  set.seed(123)
  library(qvalue)
  library(dplyr)
  trait<-"Trait.1"
  data<- .testData()

  pheno<- data$phenotype
  gene_count<- data$gene_count
  covars<-"Age,Sex"

  res<-lm_count_mat(count_matrix=gene_count,pheno=pheno,trait=trait,
               covariates_string=covars, log_transform="log_replace_half_min")

  expect_equal(colnames(res),c("geneID", "beta", "se", "t_stat", "p_value", "fdr_bh", "z_score"))
  expect_equal(dim(res),c(1000L, 7L) )

  expect_error(lm_count_mat(count_matrix=gene_count,pheno=pheno,trait=trait,
                            covariates_string=covars, log_transform="log_replace_half"),"Requested transformation not allowed.
         Allowed transformation names are log_replace_half_min, log_add_min, log_add_0.5")


  expect_identical(res$geneID, rownames(gene_count))

  expect_true(all(rownames(pheno)%in%colnames(gene_count)))
  set.seed(10)
  genes_id<-sample(rownames(gene_count),10)

  expect_message(lm_count_mat(count_matrix=gene_count,pheno=pheno,trait=trait,
                             covariates_string=covars, log_transform="log_replace_half_min", gene_IDs = genes_id),
                 "Filtering count_matrix to genes : gene.491 gene.649 gene.330 gene.368 gene.460 gene.439 gene.584 gene.438 gene.423 gene.511")

  expect_equal(nrow(lm_count_mat(count_matrix=gene_count,pheno=pheno,trait=trait,
                              covariates_string=covars, log_transform="log_replace_half_min", gene_IDs = genes_id)),length(genes_id))


})


test_that("Wrapper fast Linear regression for single exposure ", {
  set.seed(123)
  library(qvalue)
  library(dplyr)
  trait<-"Trait.1"
  data<- .testData()

  pheno<- data$phenotype
  gene_count<- data$gene_count
  covars<-"Age,Sex"

  gene_count_norm<- log_transform_count(count_matrix =gene_count, transform = "log_replace_half_min")

  filter_gene<- apply_filters(count_matrix=gene_count_norm, median_min = 1, expression_sum_min = 10,
                              max_min = 5, range_min = 5, prop_zero_max = 0.7,
                              cv_min = NULL, cv_max = NULL,
                              max_to_median_max = NULL)

  res<-lm_count_mat_emp_pval(count_matrix=filter_gene, pheno=pheno, trait=trait, covariates_string=covars,
                             n_permute=100,
                             gene_IDs=NULL,
                             log_transform = "log_replace_half_min",
                             seed = NULL,
                             stat_type="z_score",
                             empirical_type = "storey",
                             t_df = NULL,
                             family="gaussian")

  expect_equal(nrow(res), nrow(filter_gene))

  expect_identical(colnames(res),c("geneID", "beta", "se", "t_stat", "p_value", "fdr_bh", "z_score",
                                   "emp_pvals", "bh_emp_pvals"))


  expect_error(lm_count_mat_emp_pval(count_matrix=filter_gene, pheno=pheno, trait=trait, covariates_string=covars,
                                            n_permute=100,
                                            gene_IDs=NULL,
                                            log_transform = "log_replace_half_min",
                                            seed = NULL,
                                            stat_type="z_score",
                                            empirical_type = "storey",
                                            t_df = NULL,
                                            family="binomial"))


})







test_that("Wrapper fast Linear regression for single exposure ", {
  set.seed(123)
  library(qvalue)
  library(dplyr)
  trait<-"Trait.1"
  data<- .testData()

  pheno<- data$phenotype
  gene_count<- data$gene_count
  covars<-"Age,Sex"

  gene_count_norm<- log_transform_count(count_matrix =gene_count, transform = "log_replace_half_min")

  filter_gene<- apply_filters(count_matrix=gene_count_norm, median_min = 1, expression_sum_min = 10,
                              max_min = 5, range_min = 5, prop_zero_max = 0.7,
                              cv_min = NULL, cv_max = NULL,
                              max_to_median_max = NULL)

  res<-lm_count_mat_emp_pval(count_matrix=filter_gene, pheno=pheno, trait=trait, covariates_string=covars,
                             n_permute=100,
                             gene_IDs=NULL,
                             log_transform = "log_replace_half_min",
                             seed = NULL,
                             stat_type="z_score",
                             empirical_type = "storey",
                             t_df = NULL,
                             family="gaussian")

  expect_equal(nrow(res), nrow(filter_gene))

  expect_identical(colnames(res),c("geneID", "beta", "se", "t_stat", "p_value", "fdr_bh", "z_score",
                                   "emp_pvals", "bh_emp_pvals"))


  expect_error(lm_count_mat_emp_pval(count_matrix=filter_gene, pheno=pheno, trait=trait, covariates_string=covars,
                                     n_permute=100,
                                     gene_IDs=NULL,
                                     log_transform = "log_replace_half_min",
                                     seed = NULL,
                                     stat_type="z_score",
                                     empirical_type = "storey",
                                     t_df = NULL,
                                     family="binomial"))


})






test_that("fast Linear regression for multipe exposure ", {
  set.seed(123)
  library(qvalue)
  library(dplyr)
  traits<-"Trait.1,Trait.2"
  data<- .testData()

  pheno<- data$phenotype
  gene_count<- data$gene_count
  covars<-"Age,Sex"

  res<-mult_lm_count_mat(count_matrix=gene_count,pheno=pheno,traits=traits,
                    covariates_string=covars, log_transform="log_replace_half_min")

  dput(colnames(res))
  expect_equal(colnames(res),c("geneID", "beta.Trait.1", "beta.Trait.2", "t_stat", "p_value", "fdr_bh", "z_score"))
  expect_equal(dim(res),c(1000L, 7L) )

  expect_error(mult_lm_count_mat(count_matrix=gene_count,pheno=pheno,traits="Trait.1",
                            covariates_string=covars, log_transform="log_replace_half"),"Only found single trait, run linear regression for single trait instead")


  expect_identical(res$geneID, rownames(gene_count))

  expect_true(all(rownames(pheno)%in%colnames(gene_count)))
  set.seed(10)
  genes_id<-sample(rownames(gene_count),10)

  expect_message(mult_lm_count_mat(count_matrix=gene_count,pheno=pheno,traits=traits,
                              covariates_string=covars, log_transform="log_replace_half_min", gene_IDs = genes_id),
                 "Filtering count_matrix to genes : gene.491 gene.649 gene.330 gene.368 gene.460 gene.439 gene.584 gene.438 gene.423 gene.511")

  expect_equal(nrow(mult_lm_count_mat(count_matrix=gene_count,pheno=pheno,traits=traits,
                                 covariates_string=covars, log_transform="log_replace_half_min", gene_IDs = genes_id)),length(genes_id))


})






test_that("Wrapper fast Linear regression for multiple exposure ", {
  set.seed(123)
  library(qvalue)
  library(dplyr)
  traits<-"Trait.1,Trait.2"
  data<- .testData()

  pheno<- data$phenotype
  gene_count<- data$gene_count
  covars<-"Age,Sex"

  gene_count_norm<- log_transform_count(count_matrix =gene_count, transform = "log_replace_half_min")

  filter_gene<- apply_filters(count_matrix=gene_count_norm, median_min = 1, expression_sum_min = 10,
                              max_min = 5, range_min = 5, prop_zero_max = 0.7,
                              cv_min = NULL, cv_max = NULL,
                              max_to_median_max = NULL)

  res<-lm_mult_count_mat_emp_pval(count_matrix=filter_gene, pheno=pheno, traits =traits, covariates_string=covars,
                             n_permute=100,
                             gene_IDs=NULL,
                             log_transform = "log_replace_half_min",
                             seed = NULL,
                             stat_type="z_score",
                             empirical_type = "storey",
                             t_df = NULL,
                             family="gaussian")

  expect_equal(nrow(res), nrow(filter_gene))

  expect_equal(colnames(res),c("geneID", "beta.Trait.1", "beta.Trait.2", "t_stat", "p_value",
                               "fdr_bh", "z_score", "emp_pvals", "bh_emp_pvals"))
  expect_equal(dim(res),c(428L, 9L) )


  expect_error(lm_mult_count_mat_emp_pval(count_matrix=filter_gene, pheno=pheno, traits=traits, covariates_string=covars,
                                     n_permute=100,
                                     gene_IDs=NULL,
                                     log_transform = "log_replace_half_min",
                                     seed = NULL,
                                     stat_type="z_score",
                                     empirical_type = "storey",
                                     t_df = NULL,
                                     family="binomial"))


})




test_that("fast Linear regression for permutation p-value ", {
  set.seed(123)
  library(qvalue)
  library(dplyr)
  trait<-"Trait.1"
  data<- .testData()

  pheno<- data$phenotype
  gene_count<- data$gene_count
  covars<-"Age,Sex"
  set.seed(1)

  gene_count_norm<- log_transform_count(count_matrix =gene_count, transform = "log_replace_half_min")

  filter_gene<- apply_filters(count_matrix=gene_count_norm, median_min = 1, expression_sum_min = 10,
                              max_min = 5, range_min = 5, prop_zero_max = 0.7,
                              cv_min = NULL, cv_max = NULL,
                              max_to_median_max = NULL)

  gene_name<-sample(rownames(filter_gene),1)



  single_transcript<- filter_gene[gene_name,]

  residual_perm<- sapply(seq_len(1000), function(x){
    permute_resids_trait(pheno = pheno,
                         trait = trait,
                         covariates_string = covars, family = "gaussian")
  })
  perm_res<-  lm_count_mat_permute(residual_permutation=residual_perm, covariates_string=covars, pheno, single_transcript)

  expect_equal(nrow(perm_res), 1000)

  expect_true(perm_res$permute[1]=="perm_1")

  expect_identical(colnames(perm_res),c("permute", "beta", "se_beta", "test_stat", "pvalue"))

})







test_that("Wrapper fast Linear regression for permutation p-value ", {
  set.seed(123)
  library(qvalue)
  library(dplyr)
  trait<-"Trait.1"
  data<- .testData()

  pheno<- data$phenotype
  gene_count<- data$gene_count
  covars<-"Age,Sex"
  set.seed(1)

  gene_count_norm<- log_transform_count(count_matrix =gene_count, transform = "log_replace_half_min")

  filter_gene<- apply_filters(count_matrix=gene_count_norm, median_min = 1, expression_sum_min = 10,
                              max_min = 5, range_min = 5, prop_zero_max = 0.7,
                              cv_min = NULL, cv_max = NULL,
                              max_to_median_max = NULL)
  set.seed(1)
  gene_name<-sample(rownames(filter_gene),10)


  perm_res<-lm_count_mat_perm_pval(count_matrix=filter_gene, pheno=pheno, trait=trait, covariates_string=covars,
                         n_permute=1000,
                         gene_IDs=gene_name,
                         log_transform = "log_replace_half_min",
                         seed = 12,
                         family="gaussian")

  expect_equal(nrow(perm_res),length(gene_name))

  expect_equal(colnames(perm_res),c("geneID", "beta", "se", "t_stat", "p_value", "fdr_bh", "z_score",
                "perm_pval"))

  expect_identical(perm_res$geneID, gene_name)

  expect_message(lm_count_mat_perm_pval(count_matrix=filter_gene, pheno=pheno, trait=trait, covariates_string=covars,
                n_permute=10,
                gene_IDs=NULL,
                log_transform = "log_replace_half_min",
                seed = 12,
                family="gaussian"), "No list gene ID/s are found. It will run permutations for all the genes and it will take long times")


})










