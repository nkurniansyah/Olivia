
test_that("Log replace half min", {
  set.seed(0)
  data<-.testData()
  gene_count<- data$gene_count
  
  gene_count_norm<- median_normalization(gene_count)
  
  
  indx_0<-which(gene_count_norm ==0, arr.ind = T)
  indx_0<-  data.frame(indx_0)

  gene_select<- gene_count_norm[rownames(gene_count_norm)==rownames(indx_0)[1],]
  gene_select[gene_select==0]<-min(gene_select[gene_select>0])/2
  test_logval<-log2(gene_select)


  log_gene_count<- log_replace_half_min(gene_count_norm)

  expect_equal(log_gene_count[rownames(indx_0)[1],],test_logval)
  

  
})



test_that("Log add half min", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count
  
  gene_count_norm<- median_normalization(gene_count)
  
  indx_0<-which(gene_count_norm ==0, arr.ind = T)
  indx_0<-  data.frame(indx_0)
  
  gene_select<- gene_count_norm[rownames(gene_count_norm)==rownames(indx_0)[1],]
  test_logval<-log2(gene_select+min(gene_select[gene_select>0])/2)
  
  
  log_gene_count<- log_add_min(gene_count_norm)
  
  expect_equal(log_gene_count[rownames(indx_0)[1],],test_logval)
  
})





test_that("Log add half", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count
  
  gene_count_norm<- median_normalization(gene_count)
  
  indx_0<-which(gene_count_norm ==0, arr.ind = T)
  indx_0<-  data.frame(indx_0)
  
  gene_select<- gene_count_norm[rownames(gene_count_norm)==rownames(indx_0)[1],]
  test_logval<-log2(gene_select+0.5)
  
  
  log_gene_count<- log_add_0.5(gene_count_norm)
  
  expect_equal(log_gene_count[rownames(indx_0)[1],],test_logval)
  
})


test_that("log Transform Count", {
  set.seed(123)
  data<-.testData()
  gene_count<- data$gene_count
  
  gene_count_norm<- median_normalization(gene_count)
  
  log_rep_halfmin<-log_replace_half_min(gene_count_norm)
  
  log_gene_count<- log_transform_count(count_matrix =gene_count_norm, transform = "log_replace_half_min")
  expect_equal(dim(gene_count_norm),c(1000L, 20L))
  
  expect_equal(log_gene_count,log_rep_halfmin)
  
  expect_error(log_transform_count(count_matrix =gene_count, transform = "median"),"Requested transformation not allowed.
         Allowed transformation names are log_replace_half_min, log_add_min, log_add_0.5")

  
})
