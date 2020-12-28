

filter_by_genes <- function(count_matrix, gene_IDs){
  
  if (!all(gene_IDs %in% rownames(count_matrix))) {
    genes_not_avail <- gene_IDs[which(!gene_IDs %in% rownames(count_matrix))]
    message(paste("Some requested gene_IDs not in count_matrix, proceeding without using ", 
                  paste(genes_not_avail, collapse = ", ")))
  }
  genes_avail <- gene_IDs[which(gene_IDs %in% rownames(count_matrix))]
  if (length(genes_avail) == 0){
    stop("None of the requested genes are in the count_matrix, stopping...")
  }
  message(paste0("Filtering count_matrix to genes", genes_avail))
  
  count_matrix<-count_matrix[rownames(count_matrix) %in% gene_IDs, , drop = FALSE]
  return(count_matrix)
  
}