[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4398)
[![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/kbchoi/scrate)](https://hub.docker.com/r/kbchoi/scrate)

# scRATE

scRATE is a model selection-based scRNA-seq quantitation algorithm that controls overfitting and unwarranted imputation of technical zeros. We first fit UMI counts with Poisson (P), Negative-Binomial (NB), Zero-inflated Poisson (ZIP), and Zero-inflated Negative Binomial models (ZINB). These models are applicable to UMI counts directly and therefore no need of arbitrary preprocessing of counts, e.g., normalization (by scaling to, for example, CPM) or log-transformation (with pseudocounts). Our model comparison enables us to compute denoised rates of gene expression using the best model which each gene data is conforming to.


* Free software: GPLv3 license


## Installation

Installation of scRATE is simple, although it may take a while as it has to compile rstan, the related packages, and all their dependencies.

```r
> devtools::install_github('churchill-lab/scRATE')
> library(scRATE)
> version()
```
```
____ ____ ____ ____ ___ ____
[__  |    |__/ |__|  |  |___
___] |___ |  \ |  |  |  |___
                   Ver:0.1.1
```

This might hick up because of hdf5r or loomR. Then, try the following.
```r
> install.packages('hdf5r')
> devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
```

## Quick Start

Here is a quick example on how to use it. We first load `scRATE` and `loomR` package. We use loom format files to store single-cell gene expression data as it is supported in both `R` and `python`. In `scRATE`, we also provide handy features that facilitate the handling of model selection results with loom format file.

```r
> library(scRATE)
> library(loomR)
```
A data file in the example is available at <ftp://churchill-lab.jax.org/analysis/scRATE/DC-like_cells.loom>.

```r
> ds <- connect('DC-like_cells.loom')
> cntmat <- t(ds$matrix[,])
> gsymb <- ds$row.attrs$GeneID[]
> ds$close_all()
> head(gsymb)
```
```
[1] "Xkr4"    "Gm1992"  "Gm37381" "Rp1"     "Rp1.1"   "Sox17"
```

We recommend using offset (or exposure) in order to reflect the difference in depth of coverage across cells while fitting the count models. The models we compare use log link function, and therefore, the offsets should be log transformed too.

```r
> exposure <- log(colSums(cntmat))
> head(exposure)
```
```
  [1] 8.815518 8.996157 9.025816 9.037771 9.062420 9.397732
```

We will pick a gene, *Cybb* (one that we know that it fits to ZINB model significantly better than the other models), and load its UMI counts into a data frame.

```r
> gg <- 4153
> gsymb[gg]
```
```
[1] "Cybb"
```
```r
> y <- cntmat[gg,]
> y
```
```
  [1]  1  0  8  7  3  6  7  1  4  0  0  2  2 10  0  3  1  1  2  4  6  3  9  3  1
 [26]  4  1  2  3  8  0  0  1  6  2  6  3  0  0  1  0  4  1  3  0  2 17  1  2  5
 [51]  0  1  1  0  6  3  4  4  3  1  3  1  7  8  0  2  4  6  2  7  0  4  0  1  0
 [76]  2  0  8  2  0  1  4  9  3  0  2  2  3  0  3  0  0  4  1  4  3  2  3  2  2
[101]  8  0  5  8  2  3  0  3  0  3  6  0  4  2 11  5  1  2 12  2  4  2  6  8  6
[126]  3  3  1  0  2  1  1 10  2  2  0  2  2  1  0  1  0  3 11  0  7  2  2  2  3
[151]  7  0  1 10  4  1  4  5  1  0  1  3  3  0  5  4  9  4  2  0  4  1  2  4  1
[176]  2  1  4  2  1  0  1  0  1  4  6  8  9  4  3
```
```r
> gexpr <- data.frame(y=y, exposure=exposure)
```
We must have `y` and `exposure` as variables. We can also append covariates to the data frame.
```r
> gexpr <- data.frame(y=y, exposure=exposure, celltype, sex)
```

We are now ready to fit the models.
```r
> model_fit <- fit_count_models(gexpr)
```
Or we can specify our own model formula (See brms documentation) if there are covariates. Note we should not use `offset()` function in the formula.
```r
> model_fit <- fit_count_models(gexpr, 'y ~ 1 + (1|sex) + (1|celltype)')
```
```
Fitting models for Cybb
Fitting data with Poisson model...
Fitting data with Negative Binomial model...
Fitting data with Zero-Inflated Poisson model...
Compiling the C++ model

Start sampling
SAMPLING FOR MODEL 'zip' NOW (CHAIN 1).

Chain 1:
Chain 1: Gradient evaluation took 0.000177 seconds
Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.77 seconds.
Chain 1: Adjust your expectations accordingly!
Chain 1:
Chain 1:
...

Fitting data with Zero-Inflated Negative Binomial model...
Compiling the C++ model
Start sampling

SAMPLING FOR MODEL 'zinb' NOW (CHAIN 1).
Chain 1:
Chain 1: Gradient evaluation took 0.000329 seconds
Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 3.29 seconds.
Chain 1: Adjust your expectations accordingly!
Chain 1:
Chain 1:
...
```

Then we compare the models with the leave-one-out cross validation test, and select the model that best fits the data.

```r
> elpd_loo <- compare_count_models(model_fit)
> elpd_loo
```
```
       elpd_diff se_diff
model4    0.0       0.0
model2  -10.0       4.5
model3  -34.9       9.2
model1 -111.0      20.9
```
```r
> select_model(elpd_loo, margin=2)
[1] 4
```

## How to cite

K. Choi, Y. Chen, D.A. Skelly, G.A. Churchill. “Bayesian model selection reveals biological origins of zero inflation in single-cell transcriptomics.” bioRxiv. doi: <https://doi.org/10.1101/2020.03.03.974808> (2020)
