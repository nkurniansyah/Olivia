#' Quantile empirical P-values
#'
#' @param statistics A vector of calculated p-value.
#' @param null_statistics A vactor of permuted test pvalue
#'
#' @return A vector of empirical pvalues
#' @examples
#' stat<- runif(1000,0,1)
#' null_stat<- runif(100000,0,1)
#' compute_quantile_empirical_pvalues(statistics=stat,null_statistics=null_stat)
#' @export

compute_quantile_empirical_pvalues <- function(statistics, null_statistics){

  Fn <- ecdf(null_statistics)
  emp_pvalues <- Fn(statistics)

  inds <- which(emp_pvalues == 0)
  if (length(inds) > 0){
    emp_pvalues[inds] <- 1/(2*length(null_statistics))
  }

  emp_pvalues
}



#' Permutation P-values
#'
#' @param pvalue A calculated p-value from single transcript/gene
#' @param null_pval A vactor of permuted test pvalue from single transcript/gene
#'
#' @return A vector of permutation p-value
#' @examples
#' pval<- 1e-05
#' null_pval<-runif(100000,0,1)
#' permutataion_pvalues(pvalue=pval,null_pval=null_pval)
#' @export
#'


permutation_pvalues<-function(pvalue, null_pval){
  perm_pval<-sum(null_pval < pvalue)/length(null_pval)
  perm_pval
}
