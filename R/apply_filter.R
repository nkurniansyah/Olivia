#'Filter gene counts matrix by applying various filters
#'
#'The purpose of filtering is to remove lowly express genes
#'
#' @param count_matrix  A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#' @param cv_min  minimum allowed value of coefficient of variation. Default is NULL.
#' @param cv_max  maximum allowed value of coefficient of variation. Default is NULL.
#' @param expression_sum_min  minimum allowed sum of counts across individuals. Default is 10.
#' @param median_min  minimum allowed median count value in the sample. Default is 1.
#' @param max_to_median_max maximum allowed ratio between the maximum to the
#' @param max_min minimum allowed maximum value in the sample for a gene. Default is 5.
#' @param range_min minimum allowed range of counts in the sample. Default is 5.
#' @param prop_zero_max maximum allowed proportion of individuals with count of
#'                      zero for a transcript. Default is 0.8.
#'
#' @return Matrix of gene counts, with filtered (reduced) rows.
#' @examples
#'
#' data(rnaseq_count_matrix)
#' 
#' apply_filters(count_matrix=rnaseq_count_matrix, median_min = 1, expression_sum_min = 10,max_min = 5,
#'               range_min = 5, prop_zero_max = 0.8,cv_min = NULL, cv_max = NULL,
#'               max_to_median_max = NULL)
#'
#' @export

apply_filters <- function(count_matrix, median_min = 1, expression_sum_min = 10,
                          max_min = 5, range_min = 5, prop_zero_max = 0.8,
                          cv_min = NULL, cv_max = NULL,
                          max_to_median_max = NULL){

  message(paste("applying filters on a transcript count matrix of",
                nrow(count_matrix), "transcripts, across", ncol(count_matrix), "individuals"))

  # check missingness
  na_inds <- which(is.na(count_matrix))
  if (length(na_inds) > 0){
    messsage(paste("There are", length(na_inds), "missing values, setting them to zero"))
    count_matrix[na_inds] <- 0
  }


  ############ Compute characteristics of the normalized count matrix ##############

  message("Computing transtripts characteristics...")
  median_vals <- apply(count_matrix, 1, median)
  expression_sum_vals <- rowSums(count_matrix)
  max_vals <- apply(count_matrix,1,max)
  range_vals <- apply(count_matrix, 1, function(x) max(x) - min(x))
  prop_zero_vals <- apply(count_matrix, 1, function(x) mean(x == 0))
  max_to_median_vals <- apply(count_matrix, 1, function(x) max(x)/median(x))
  cv_vals <- apply(count_matrix, 1, function(x) sd(x)/mean(x))

  message("Appying filters...")
  inds_rm <- c()
  if (!is.null(median_min)){
    inds_median_min <- which(median_vals < median_min)
    message(paste("There are", length(inds_median_min), "transcripts with median
                  value lower than", median_min))
    inds_rm <- c(inds_rm, inds_median_min)
  }

  if (!is.null(expression_sum_min)){
    inds_expression_sum_min <- which(expression_sum_vals < expression_sum_min)
    message(paste("There are", length(inds_expression_sum_min), "transcripts with expression sum
                  value lower than", expression_sum_min))
    inds_rm <- c(inds_rm, inds_expression_sum_min)
  }

  if (!is.null(max_min)){
    inds_max_min <- which(max_vals < max_min)
    message(paste("There are", length(inds_max_min), "transcripts with maximum expression
                  value lower than", max_min))
    inds_rm <- c(inds_rm, inds_max_min)
  }

  if (!is.null(range_min)){
    inds_range_min <- which(range_vals < range_min)
    message(paste("There are", length(inds_range_min), "transcripts with maximum
                  expression range value lower than", range_min))
    inds_rm <- c(inds_rm, inds_range_min)
  }

  if (!is.null(prop_zero_max)){
    inds_prop_zero_max <- which(prop_zero_vals > prop_zero_max)
    message(paste("There are", length(inds_prop_zero_max), "transcripts with propotion
                  of zero counts higher than", prop_zero_max))
    inds_rm <- c(inds_rm, inds_prop_zero_max)
  }

  if (!is.null(cv_min)){
    inds_cv_min <- which(cv_vals < cv_min)
    message(paste("There are", length(inds_cv_min), "transcripts with coefficient of
                  variation lower than", cv_min))
    inds_rm <- c(inds_rm, inds_cv_min)
  }

  if (!is.null(cv_max)){
    inds_cv_max <- which(cv_vals < cv_max)
    message(paste("There are", length(inds_cv_max), "transcripts with coefficient of
                  variation higher than", cv_max))
    inds_rm <- c(inds_rm, inds_cv_max)
  }

  if (!is.null(max_to_median_max)){
    inds_max_to_median_max <- which(max_to_median_vals > max_to_median_max)
    message(paste("There are", length(inds_max_to_median_max), "transcripts with max/median  value
                   higher than", max_to_median_max))
    inds_rm <- c(inds_rm, inds_max_to_median_max)
  }

  inds_rm <- unique(inds_rm)
  message(paste("Removing", length(inds_rm), "unique transcripts not passing requested filters"))


  count_matrix <- count_matrix[-inds_rm, ]
  return(count_matrix)
}
