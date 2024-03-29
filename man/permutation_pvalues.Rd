% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_empirical_pvalue.R
\name{permutation_pvalues}
\alias{permutation_pvalues}
\title{Permutation P-values}
\usage{
permutation_pvalues(pvalue, null_pval)
}
\arguments{
\item{pvalue}{A calculated p-value from single transcript/gene}

\item{null_pval}{A vactor of permuted test pvalue from single transcript/gene}
}
\value{
A vector of permutation p-value
}
\description{
Permutation P-values
}
\examples{
pval<- 1e-05
null_pval<-runif(100000,0,1)
permutataion_pvalues(pvalue=pval,null_pval=null_pval)
}
