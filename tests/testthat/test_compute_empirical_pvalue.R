test_that("Quantile empirical p-values", {
  set.seed(123)
  stat<- runif(1000,0,1)
  null_stat<- runif(100000,0,1)
  emp_pval<-compute_quantile_empirical_pvalues(statistics=stat,
                                               null_statistics=null_stat)
  expect_equal(length(stat),length(emp_pval))
})


# we don't need it-- we haven't developed it, we are just providing a wrapper.
test_that("Storey empirical p-values", {
  set.seed(123)
  library(qvalue)
  stat<- rnorm(1000,0,1)
  null_stat<- rnorm(100000,0,1)
  sto_emp_pval<-compute_storey_empirical_pvalues(statistics=stat,null_statistics=null_stat)
  expect_equal(length(stat),length(sto_emp_pval))
})



test_that("Permutation p-values", {
  pval<- 1e-05
  null_pval<-runif(100000,0,1)
  perm_pval<-permutation_pvalues(pvalue=pval,null_pval=null_pval)
  expect_equal(length(pval),length(perm_pval))
  expect_true( abs(pval - perm_pval) < 1e-4)
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
  
  expect_equal(empval,compute_quantile_empirical_pvalues(statistics=stat,
                                                         null_statistics=null_stat) )
  
})
