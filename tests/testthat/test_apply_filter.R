test_that("Apply filters", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count
  
  gene_count_norm<- median_normalization(gene_count)
  
  # use the apply_filters command for each characteristic separately. 
  
  filtered_gene<- apply_filters(count_matrix=gene_count_norm, median_min = 1, expression_sum_min = 10,
                                max_min = 5, range_min = 5, prop_zero_max = 0.5,
                                cv_min = 0, cv_max = 0.5,
                                max_to_median_max = 4)
  
  # we want to see that after filtering, the data satisfies certain properties.
  
  # make sure median number of transcripts is always at least 1: 
  expect_true(all(apply(filtered_gene, 1, median) > 0))
  
  expect_true(all(rowSums(filtered_gene)>10))
  expect_true(all(apply(filtered_gene, 1, max) > 5))
  expect_true(all(apply(filtered_gene, 1, function(x) max(x) - min(x)) > 5))
  
  expect_true(all(apply(filtered_gene, 1, function(x) mean(x == 0)) < 0.5))
  
  expect_true(all(apply(filtered_gene, 1,function(x) max(x)/median(x)) < 4))

  expect_true(all(apply(filtered_gene, 1, function(x) sd(x)/mean(x))>0.5))
  
  
})