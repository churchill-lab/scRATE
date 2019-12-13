======
scRATE
======

scRATE is a model selection-based scRNA-seq quantitation algorithm that controls overfitting and unwarranted imputation of technical zeros. We first fit UMI counts with Poisson, Negative-Binomial, Zero-inflated Poisson, and Zero-inflated Negative Binomial models. These models are applicable to UMI counts directly and therefore no need of arbitrary preprocessing of counts, e.g., normalization (by scaling to, for example, CPM) or log-transformation (with pseudocounts). Our model comparison enables us to compute denoised rates of gene expression using the best model which each gene data is conforming to.


* Free software: MIT license
* Documentation: under development


How to cite
-----------

* Manuscript in preparation (expected in early 2020)

