test_that("Quantile empirical p-values", {
  set.seed(123)
  stat<- runif(1000,0,1)
  null_stat<- runif(100000,0,1)
  emp_pval<-compute_quantile_empirical_pvalues(statistics=stat,
                                               null_statistics=null_stat)
  expect_equal(length(stat),length(emp_pval))
})





test_that("Permutation p-values", {
  pval<- 1e-05
  null_pval<-runif(100000,0,1)
  perm_pval<-permutation_pvalues(pvalue=pval,null_pval=null_pval)
  expect_equal(length(pval),length(perm_pval))
  expect_true( abs(pval - perm_pval) < 1e-4)
})


