#' Filter gene counts matrix by geneID
#'
#' Filtering gene counts matrix based on intrested gene/s
#'
#' @param count_matrix  A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#' @param gene_IDs A vactor of gene names
#'
#' @return Matrix of selected gene counts
#' @examples
#' data(rnaseq_count_matrix)
#' genes<-c("ENSG00000000003","ENSG00000000005","ENSG00000000419")
#' filter_by_genes(count_matrix=rnaseq_count_matrix, gene_IDs=genes)
#' @export

filter_by_genes <- function(count_matrix, gene_IDs){

  if (!all(gene_IDs %in% rownames(count_matrix))) {
    genes_not_avail <- gene_IDs[which(!gene_IDs %in% rownames(count_matrix))]
    message(paste("Some requested gene_IDs not in count_matrix, proceeding without using ",
                  paste(genes_not_avail, collapse = ", ")))
  }
  genes_avail <- gene_IDs[which(gene_IDs %in% rownames(count_matrix))]
  if (length(genes_avail) == 0){
    stop("None of the requested genes are in the count_matrix, stopping..")
  }

  message(paste(c("Filtering count_matrix to genes :", genes_avail), collapse= " "))

  count_matrix<-count_matrix[rownames(count_matrix) %in% gene_IDs, , drop = FALSE]
  return(count_matrix)
}
