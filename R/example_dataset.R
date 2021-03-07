#' Simulated phenotypes.
#'
#' A dataset containing simultaed 40 samples
#'
#' @format A data frame with 40 rows and 4 variables, and row names must be matched with column names in rnaseq_count_matrix:
#' \describe{
#'   \item{Age}{covariate, must be integer}
#'   \item{Sex}{covariate,0 is female and 1 is male, binary}
#'   \item{Race}{covariate,race must be categorical}
#'   \item{Trait.1}{outcome/ trait , numeric}
#'   \item{Trait.2}{outcome/ trait , numeric}
#'   ...
#' }

"phenotype"


#' Raw gene count matrix
#'
#' Example of bulk RNA-seq
#'
#' @format a matrix which contains 58051 and 40 column.
#' \describe{
#'   \item{column}{sample names and it must be matched with the rownames in the phenotype}
#'   \item{rownames}{ensemblID/ geneID }

#'  @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151243}
#' }

"rnaseq_count_matrix"




#' Single gene count (ENSG00000000003)
#'
#' Example of single gene count
#'
#' @format a vector which contains  40 column and match with rownames in phenotypes
#' \describe{
#'   \item{column}{sample names and it must be matched with the rownames in the phenotype}

#'  @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151243}
#' }

"ENSG00000000003"
