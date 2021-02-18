
test_that("Median Normalization", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count
  median_expression_sum <- median(colSums(gene_count))
  
  med_norm<- median_normalization(gene_count)
  expect_true(is.numeric(med_norm), TRUE)
  expect_equal(all(colSums(med_norm) == median_expression_sum), TRUE)
  
})