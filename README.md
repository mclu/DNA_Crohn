# Analysis of Methylation for Crohn's Disease

## Overview
The repo analyzes DNA methylation in Adipocyte Tissues for Crohn's Disease, along with demonstrating the usage of permutation tests and parallel computing.

## Navigation
- 
- 

## Installation
Data can be downloaded using the command line. Details of the data and experiment can be found on the [Gene Expression Omnibus website](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138311)

```bash
wget -O - ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138311/matrix/GSE138311_series_matrix.txt.gz | gunzip -c > matrix.txt
```

To run the R script, the following packages should be pre-installed in the IDE.
```r
install.packages('tidyverse')
install.packages('data.table')
install.packages('parallel')
install.packages('future')
```

## Remark
The repo is  revised work of the course [Stats506 Computational Methods and Tools in Statistics](https://jbhender.github.io/Stats506/F19/index.html) taught by Dr. Henderson.
