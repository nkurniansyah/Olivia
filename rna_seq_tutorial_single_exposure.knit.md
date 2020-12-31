---
title: "Benchmarking Association Analyses of Continuous Exposures with RNA-seq in Observational Studies: tutorial"
author: "Tamar Sofer & Nuzulul Kurniansyah"
date: "12/29/2020"
output: 
  html_document:
    toc: true
  pdf_document:
    toc: true
---



# Introduction
Here we demonstrate how to perform association analyses of continuous phenotypes 
with RNA-seq data based on the pipeline proposed in the manuscript ``Benchmarking Association Analyses of Continuous Exposures with RNA-seq in Observational Studies". 


```r
base_path <- "~/Documents/GitHub/RNASeq/"
```

# Setting up: packages, functions

## Installation Requirements

- An installed R distribution, such as Microsoft R Open. Use R version 3.6.3 or later.
- RStudio, either the commercial edition or the open-source RStudio Desktop.
- R packages (CRAN)



```r
list.of.packages <- c("dplyr", "tidyverse","purrr","data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dplyr)
library(tidyverse)
library(purrr)
library(data.table)
```

Install packages from Bioconductor
(we may remove some of these, if we are not doing any comparisons in this code)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

list.of.packages <- c("limma", "qvalue","DESeq2","edgeR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
library(limma)
library(qvalue)
library(DESeq2)
library(edgeR)
```

Source files with functions (this will change once the code is packaged)

```r
files <- list.files(paste0(base_path, "Code"), full.names = TRUE)
for (file in files)(
  source(file)
)
```


# Load example data
## Load transcripts 
First we load the transcripts, where were obtained from XXXX.
We reformat to a desired form. 


```r
rnaseq_matrix <- fread(paste0(base_path, "Data/GSE156437_RNA_seq.txt"), data.table = F)

# make transcript names rownames rather than a column
row.names(rnaseq_matrix) <- rnaseq_matrix[,"ensg"]
rnaseq_matrix <- rnaseq_matrix[,-grep("ensg", colnames(rnaseq_matrix))]

# make a matrix instead of a data.frame
rnaseq_matrix <- as.matrix(rnaseq_matrix)

# take a look 
dim(rnaseq_matrix)
```

```
## [1] 50692    46
```

```r
rnaseq_matrix[1:10,1:3]
```

```
##                    OKS_shLuc_D1_re1 OKS_shLuc_D1_re2
## ENSMUSG00000000001           110.14            71.53
## ENSMUSG00000000003             0.00             0.00
## ENSMUSG00000000028            19.66            31.78
## ENSMUSG00000000031           150.21           119.75
## ENSMUSG00000000037             1.86             0.62
## ENSMUSG00000000049             0.00             0.00
## ENSMUSG00000000056             5.77             8.69
## ENSMUSG00000000058            50.51            48.84
## ENSMUSG00000000078            73.13            50.71
## ENSMUSG00000000085            19.94            47.73
##                    OKS_shYthdf2_D1_re1
## ENSMUSG00000000001              106.90
## ENSMUSG00000000003                0.00
## ENSMUSG00000000028               10.62
## ENSMUSG00000000031              303.06
## ENSMUSG00000000037                1.51
## ENSMUSG00000000049                0.00
## ENSMUSG00000000056                7.82
## ENSMUSG00000000058               52.76
## ENSMUSG00000000078               76.55
## ENSMUSG00000000085               46.36
```
After looking at a few transcript counts, we infer that the data was normalized already, because the counts are not natural numbers. We ignore this fact, force the entries to be integerse, and continue demonstrating the pipeline. 


```r
rnaseq_matrix <- round(rnaseq_matrix)
storage.mode(rnaseq_matrix) <- "integer"
```

## Load simulated phenotypes
We simulated in advance a data.frame phenotypes. 


```r
phenotypes <- fread(paste0(base_path, "Data/2020-08-28_pheno.csv"), data.table = F)
head(phenotypes)
```

```
##                    ID Age      BMI Sex     Trait.1
## 1    OKS_shLuc_D1_re1  35 15.83645   1   3.3999306
## 2    OKS_shLuc_D1_re2  45 54.09436   0 111.1469390
## 3 OKS_shYthdf2_D1_re1  77 17.50906   0  -1.0433947
## 4 OKS_shYthdf2_D1_re2  57 82.06248   0  87.8721591
## 5 OKS_shYthdf3_D1_re1  60 17.07899   0  -0.4006099
## 6 OKS_shYthdf3_D1_re2  32 49.14698   0 103.2596224
##     Trait.2  Race
## 1  14.01064 Black
## 2  87.81247 White
## 3  12.21969 Black
## 4 176.66123 White
## 5  16.45065 Black
## 6  55.44270 Asian
```

We define the trait of interest to study as an exposure associated with transcription. It has to be a column name in the phenotype data.frame.


```r
trait <- "Trait.1"
```

We will adjust our analysis to the simulated covariates Age, BMI, and Sex. The covriates have to correspond to column names in the phenotype data.frame. In the analysis, we will use a string defining the regression model (just the covariates part of it), so we define it here: 


```r
covariates_string <- "Age + BMI + Sex"
```
Note that we can also define the string to be "Age + BMI + as.factor(Sex)", or use interaction terms, like one would use in regression functions in R. 

Match the (simulated) individuals between the phenotype and the RNA-seq count matrix. Make sure the IDs overlap.

```r
IDs_both <- intersect(phenotypes$ID, colnames(rnaseq_matrix))
rnaseq_matrix <- rnaseq_matrix[, IDs_both] 
phenotypes <- phenotypes[match(IDs_both, phenotypes$ID),]
```
 

# Normalize the RNA-seq dataset
The options of normalization method are: median_normalization, SizeFactor, TMM. 
If the selected normalization method is SizeFactor, the users needs to define the phenotype, outcome, and covariate string. Here we show how each of these is used. We move forward in this tutorial with the median normalization, because it was used in the manuscript. However, there are no downstream differences in how the methods are applied once the data is normalized. 



```r
median_norm <- normalize_trancript_count(rnaseq_matrix,
                                          normalization_type = "median_normalization")
```

```
## Applying median normalization...
```

```r
tmm_norm <- normalize_trancript_count(rnaseq_matrix,
                                      normalization_type = "TMM")
```

```
## Applying TMM normalization...
```

```r
sizefactor_norm <- normalize_trancript_count(rnaseq_matrix,
                                             phenotypes=phenotypes, 
                                             covariates_string=covariates_string, 
                                             trait=trait, 
                                             normalization_type = "SizeFactor")
```

```
## Applying size factor normalization...
```

```
##   the design formula contains one or more numeric variables with integer values,
##   specifying a model with increasing fold change for higher values.
##   did you mean for this to be a factor? if so, first convert
##   this variable to a factor using the factor() function
```

```
##   the design formula contains one or more numeric variables that have mean or
##   standard deviation larger than 5 (an arbitrary threshold to trigger this message).
##   it is generally a good idea to center and scale numeric variables in the design
##   to improve GLM convergence.
```


## Filtering transcripts

Remove lowly express transcripts.



```r
clean_count_matrix <- apply_filters(count_matrix = median_norm, 
                                   median_min = 1, 
                                   expression_sum_min = 10, 
                                   max_min = 10, 
                                   range_min = 5, 
                                   prop_zero_max = 0.5)
```

```
## applying filters on a transcript count matrix of 50692 transcripts, across 46 individuals
```

```
## Computing transtripts characteristics...
```

```
## Appying filters...
```

```
## There are 36855 transcripts with median 
##                   value lower than 1
```

```
## There are 29223 transcripts with expression sum 
##                   value lower than 10
```

```
## There are 35758 transcripts with maximum expression
##                   value lower than 10
```

```
## There are 32769 transcripts with maximum 
##                   expression range value lower than 5
```

```
## There are 35777 transcripts with propotion 
##                   of zero counts higher than 0.5
```

```
## Removing 38236 unique transcripts not passing requested filters
```


After filtering transcripts, there are 12456 remaining for differential expression analysis. 


# Perform differential expression analysis
First, we show how we perform differential expression analysis on all transcripts without computing empirical p-values. 

```r
deg <- lm_count_mat(count_matrix=clean_count_matrix, 
                    pheno=phenotypes, 
                    trait=trait, 
                    covariates_string=covariates_string,
                    gene_IDs=NULL, 
                    log_transform = "log_replace_half_min")
head(deg)
```

```
##               geneID          beta          se
## 1 ENSMUSG00000000001 -0.0045504760 0.001739805
## 2 ENSMUSG00000000028  0.0016100430 0.002531851
## 3 ENSMUSG00000000031 -0.0074099288 0.012456586
## 4 ENSMUSG00000000056 -0.0022535102 0.002800340
## 5 ENSMUSG00000000058 -0.0057464967 0.006161747
## 6 ENSMUSG00000000078 -0.0004788995 0.004638027
##       t_stat    p_value    FDR_BH    z_score
## 1 -2.6155085 0.01241258 0.3280221 -2.5001926
## 2  0.6359154 0.52836557 0.9415882  0.6305030
## 3 -0.5948603 0.55520584 0.9482527 -0.5899773
## 4 -0.8047274 0.42561902 0.9112031 -0.7967108
## 5 -0.9326083 0.35648204 0.8912767 -0.9220892
## 6 -0.1032550 0.91826400 0.9947777 -0.1026207
```
Let's look at the distribution of p-values that we got: 

```r
hist(deg$p_value)
```

<img src="rna_seq_tutorial_single_exposure_files/figure-html/unnamed-chunk-14-1.png" width="672" />


## Use residual permutation to compute p-values under the null
We want to generate p-values under the null, in order to compute empirical p-values. Therefore, we create a ``residual permuted" trait 100 times, perform differential expression analysis, and use the resulting p-values as our null p-values. 



```r
n_permute <- 100
permuted_trait <- matrix(NA, nrow = nrow(phenotypes), ncol = n_permute)

# generated permuted traits (by residual permutation)
set.seed(12)
for (i in 1:n_permute){
  permuted_trait[,i] <- permute_resids_trait(pheno = phenotypes,
                                                   trait = trait,
                                                   covariates_string = covariates_string)
}

# perform differential expression analysis for each of the permuted traits and save z-scores
# and p-values (we don't really need both, this is for demonstration)
null_z_scores <- null_pvals <- matrix(NA, nrow = nrow(deg), ncol = n_permute)

for (i in 1:n_permute){
  phenotypes$perm_trait <- permuted_trait[,i]
  deg_perm <- lm_count_mat(count_matrix=clean_count_matrix, 
                    pheno=phenotypes, 
                    trait="perm_trait", 
                    covariates_string=covariates_string,
                    gene_IDs=NULL, 
                    log_transform = "log_replace_half_min")
  null_z_scores[,i] <- deg_perm$z_score
  null_pvals[,i] <- deg_perm$p_value
}
```

# Computing emprical p-values
We demonstrate computing quantile empirical p-values and Storey empirical p-values (see manuscript).


```r
quantile_empirical_pvals <- compute_empirical_pvalues(
                                  statistics = deg$p_value,
                                  null_statistics = as.numeric(null_pvals),
                                  stat_type = "p_value", 
                                  empirical_type = "quantile")

storey_empirical_pvals <- compute_empirical_pvalues(
                                  statistics = deg$z_score,
                                  null_statistics = as.numeric(null_z_scores),
                                  stat_type = "z_score",
                                  empirical_type = "storey"
)
```

The distribution of the empirical p-values based on the quantile method:

```r
hist(quantile_empirical_pvals)
```

<img src="rna_seq_tutorial_single_exposure_files/figure-html/unnamed-chunk-17-1.png" width="672" />

The distribution of the empirical p-values based on Storey's implementation: 

```r
hist(storey_empirical_pvals)
```

<img src="rna_seq_tutorial_single_exposure_files/figure-html/unnamed-chunk-18-1.png" width="672" />

Alternatively, this analysis can be done in a single function 
(here, stat_type is always z_score)


```r
res <- lm_count_mat_emp_pval(count_matrix = clean_count_matrix, 
                             n_permute = n_permute, 
                             seed = 12, 
                             pheno=phenotypes, 
                             trait = trait,
                             covariates_string=covariates_string,
                             gene_IDs=NULL, 
                             log_transform = "log_replace_half_min", 
                             empirical_type = "storey")
```

```
## Performing residual permutation to generate permuted trait...
```

```
## performing differential expression analysis on 100 permuted traits
```

```
## Computing empirical p-values
```
