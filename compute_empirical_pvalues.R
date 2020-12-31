
#' Title: Compute empirical P-values
#'
#' @param statistics  A vector of statistics (p-value, z-scores, or t-statistics) 
#' @param null_statistics  A vector of statistics under null, same type as statistics
#' @param empirical_type Type of requested empirical p-values. quantile or storey.
#' @param empirical_type Type of provided statistics, could be p_value, z_score, or t_stat
#' @param t_df Degrees of freedom of t-stat if provided.  
#'
#' @return vector of empirical pvalues

compute_empirical_pvalues <- function(statistics, 
                                      null_statistics, 
                                      empirical_type = "quantile",
                                      stat_type = "p_value", 
                                      t_df = NULL){
  
  if (!is.element(empirical_type, c("quantile", "storey"))){
    stop(paste("Requested empirical p-value types (empirical_type) is", empirical_type,
           "allowed values are quantile and storey", "\n"))
  }
  
  if (!is.element(stat_type, c("z_score", "p_value", "t_stat"))){
    stop(paste("Provided statistics (stat_type) are", stat_type,
               "allowed types are p_value, z_score, and t_stat", "\n"))
  }
  
  if (empirical_type == "quantile"){
    if (stat_type == "z_score"){
      statistics <- pchisq(statistics^2, df =1, lower.tail = FALSE)
      null_statistics <- pchisq(null_statistics^2, df =1, lower.tail = FALSE)
    } 
    
    if (stat_type == "t_stat"){
      if (is.null(t_df)) stop("Missing degress of freedome (t_df) for t-stat")
      statistics <- 2*pt(abs(statistics), df = t_df, lower.tail = FALSE)
      null_statistics <- 2*pt(abs(null_statistics), df = t_df, lower.tail = FALSE)
    }

    compute_quantile_empirical_pvalues(statistics, null_statistics)
  } else{  # empirical_type == "storey"
    if (stat_type == "p_value") 
        stop("Storey empirical p-value need z-scores or t-statistics")
    
    compute_storey_empirical_pvalues(statistics, null_statistics)
  }
  
  
}




#' Title: Quantile empirical P-values
#'
#' @param p_values 
#' @param null_p_values 
#'
#' @return vector of empirical pvalues

compute_quantile_empirical_pvalues <- function(p_values, null_p_values){
  
  Fn <- ecdf(null_p_values)
  emp_pvalues <- Fn(p_values)
  
  inds <- which(emp_pvalues == 0)
  if (length(inds) > 0){
    emp_pvalues[inds] <- 1/(2*length(null_p_values))
  }
  
  emp_pvalues
}



compute_storey_empirical_pvalues <- function(z_score, null_z_score){
  emp_pvalues <- empPvals(z_score, null_z_score)
  emp_pvalues
}

