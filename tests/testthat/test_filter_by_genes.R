test_that("filter by gene", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count
  
  n_sample_genes <- 5
  genes<- sample(rownames(gene_count), n_sample_genes)
  expect_equal(length(genes), n_sample_genes)
  expect_true(all(genes  %in% rownames(gene_count)))
  
  gene_selected<- filter_by_genes(count_matrix = gene_count, gene_IDs =genes )
  
  expect_equal(nrow(gene_selected),n_sample_genes)
  
  expect_equal(all(rownames(gene_selected) %in% genes), TRUE)
  
})
