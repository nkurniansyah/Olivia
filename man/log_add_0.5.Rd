% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_transform.R
\name{log_add_0.5}
\alias{log_add_0.5}
\title{Log add half}
\usage{
log_add_0.5(count_matrix)
}
\arguments{
\item{count_matrix}{A matrix of gene counts (possibly transformed). rows are genes, columns are individuals}
}
\value{
A matrix of gene expression counts after transformation
}
\description{
Type of log transform, adding 0.5 into entire gene counts matrix
}
\examples{
data(rnaseq_count_matrix)
rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
log_add_0.5(count_matrix=rnaseq_count_matrix)
}
