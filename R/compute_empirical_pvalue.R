#' Computing empirical pvalue
#'
#' @param statistics A vector of calculated test statistics
#' @param null_statistics A vactor of permuted test statistic
#' @param empirical_type Type of empirical pvalue : quantile or storey
#' @param stat_type Statistic type : p_value (quantile empirical pvalue), t_stat and chisq_stat (storey)
#' @param t_df A vector of calculated t-statistic
#' @return A vector of empirical p-value
#' @examples
#' library(qvalue)
#' stat<- runif(1000,0,1)
#' null_stat<- runif(100000,0,1)
#' compute_empirical_pvalues(statistics=stat,null_statistics=null_stat,
#'                          empirical_type = "quantile",stat_type = "p_value",
#'                          t_df = NULL)
#' @export


compute_empirical_pvalues <- function(statistics,
                                      null_statistics,
                                      empirical_type = "storey",
                                      stat_type = "t_stat",
                                      t_df = NULL){

  if (!is.element(empirical_type, c("quantile", "storey"))){
    stop(paste("Requested empirical p-value types (empirical_type) is", empirical_type,
               "allowed values are quantile and storey", "\n"))
  }

  if (!is.element(stat_type, c("chisq_stat", "p_value", "t_stat"))){
    stop(paste("Provided statistics (stat_type) are", stat_type,
               "allowed types are p_value, and t_stat", "\n"))
  }

  message(paste0("Run ",empirical_type," empirical p-values ","using ",stat_type))

  if (empirical_type == "quantile"){

    if (stat_type == "t_stat"){
      if (is.null(t_df)) stop("Missing degress of freedom (t_df) for t-stat")
      statistics <- 2*pt(abs(statistics), df = t_df, lower.tail = FALSE)
      null_statistics <- 2*pt(abs(null_statistics), df = t_df, lower.tail = FALSE)
    }
    
    compute_quantile_empirical_pvalues(statistics, null_statistics)
    
  } else{  # empirical_type == "storey"
    if (stat_type == "p_value")
      stop("Storey empirical  need chisq-statistics or t-statistics")

    compute_storey_empirical_pvalues(statistics, null_statistics)
  }


}




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


#' Storey empirical P-values
#'
#' @param statistics A vector of calculated test statistics
#' @param null_statistics A vactor of permuted test statistic
#'
#' @return A vector of empirical pvalues
#'
#' @references
#' Storey J, Bass A, Dabney A, Robinson D. 2019. qvalue: Q-value estimation for false discovery rate control. In R package version 2.18.0
#'
#' @examples
#' library(qvalue)
#' stat<- rnorm(1000,0,1)
#' null_stat<- rnorm(100000,0,1)
#' compute_storey_empirical_pvalues(statistics=stat,null_statistics=null_stat)
#'
#' @export

compute_storey_empirical_pvalues <- function(statistics, null_statistics){
  emp_pvalues <- empPvals(statistics, null_statistics)
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
