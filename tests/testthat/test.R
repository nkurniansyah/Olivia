

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

  gene_count_norm<- log_transform_count(count_matrix =gene_count, 
                                        transform = "log_replace_half_min")

  filter_gene<- apply_filters(count_matrix=gene_count_norm, 
                              median_min = 1, expression_sum_min = 10,
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










