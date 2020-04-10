# EpiTxDb <img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/EpiTxDb/EpiTxDb.png" height="200" align="right">

<!-- badges: start -->
[![R-CMD-check](https://github.com/FelixErnst/EpiTxDb/workflows/R-CMD-check/badge.svg)](https://github.com/FelixErnst/EpiTxDb/actions/)
[![Build Status](https://travis-ci.com/FelixErnst/EpiTxDb.svg?branch=master)](https://travis-ci.com/FelixErnst/EpiTxDb)
[![BioC Build](https://bioconductor.org/shields/build/devel/bioc/EpiTxDb.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/EpiTxDb/)
[![codecov](https://codecov.io/gh/FelixErnst/EpiTxDb/branch/master/graph/badge.svg)](https://codecov.io/gh/FelixErnst/EpiTxDb)
[![BioC Years](https://bioconductor.org/shields/years-in-bioc/EpiTxDb.svg)](https://doi.org/doi:10.18129/B9.bioc.EpiTxDb)
<!-- badges: end -->


The epitranscriptome includes all post-transcriptional modifications of the RNA
and describes and additional layer of information encoded on RNA. Like the term
epigenome it is not about a change in nucleotide sequences, but the addition of
functional elements through modifications.

With the development of high throughput detection strategies for specific RNA
modifications, such as miCLIP and Pseudo-Seq amongst other, a large number of
modified positions have been identified and were summarized via the RMBase 
project ([Xuan et al. 2017, Sun et al. 2015](#Literature)) project.

To make these information avaialble within the Bioconductor universe `EpiTxDb`
was developed, which facilitates the storage of epitranscriptomic information.
More specifically, it can keep track of modification identity, position, the
enzyme for introducing it on the RNA, a specifier which determines the position
on the RNA to be modified and the literature references each modification is
associated with.

# Installation

The current version of the `EpiTxDb` package is available from Bioconductor.

## Github

```
remotes::install_github("FelixErnst/EpiTxDb")
#
library(EpiTxDb)
```

## Bioconductor

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')
BiocManager::install("EpiTxDb")
library(EpiTxDb)
```

# Literature

- Jia-Jia Xuan, Wen-Ju Sun, Ke-Ren Zhou, Shun Liu, Peng-Hui Lin, Ling-Ling
Zheng, Liang-Hu Qu, Jian-Hua Yang (2017): "RMBase v2.0: Deciphering the Map of
RNA Modifications from Epitranscriptome Sequencing Data." Nucleic Acids
Research, Volume 46, Issue D1, 4 January 2018, Pages D327–D334.
doi:[10.1093/nar/gkx934](https://doi.org/10.1093/nar/gkx934)

- Wen-Ju Sun, Jun-Hao Li, Shun Liu, Jie Wu, Hui Zhou, Liang-Hu Qu, Jian-Hua Yang
(2018): "RMBase: a resource for decoding the landscape of RNA modifications from
high-throughput sequencing data", Nucleic Acids Research, Volume 44, Issue D1, 4
January 2016, Pages D259–D265.
doi:[10.1093/nar/gkv1036](https://doi.org/10.1093/nar/gkv1036).
