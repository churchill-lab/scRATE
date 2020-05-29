# scRATE

scRATE is a model selection-based scRNA-seq quantitation algorithm that controls overfitting and unwarranted imputation of technical zeros. We first fit UMI counts with Poisson (P), Negative-Binomial (NB), Zero-inflated Poisson (ZIP), and Zero-inflated Negative Binomial models (ZINB). These models are applicable to UMI counts directly and therefore no need of arbitrary preprocessing of counts, e.g., normalization (by scaling to, for example, CPM) or log-transformation (with pseudocounts). Our model comparison enables us to compute denoised rates of gene expression using the best model which each gene data is conforming to.


* Free software: GPLv3 license
* Documentation: under development


## Installation

Installation of scRATE is simple, although it may take a while as it compiles rstan and the related packages.

```r
> devtools::install_github('churchill-lab/scRATE')
> library(scRATE)
> version()

┌─┐┌─┐╦═╗╔═╗╔╦╗╔═╗
└─┐│  ╠╦╝╠═╣ ║ ║╣ 
└─┘└─┘╩╚═╩ ╩ ╩ ╚═╝
         Ver:0.1.0
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
A data file in the example is available <a href=">here</a>.

```r
> ds <- connect('DC-like_cells.loom')
> cntmat <- ds$matrix[,]
> gsymb <- ds$row.attrs$GeneID[]
> head(gsymb)
[1] "Xkr4"    "Gm1992"  "Gm37381" "Rp1"     "Rp1.1"   "Sox17"
> ds$close_all()
```

We recommend using offset (or exposure) in order to reflect the difference in depth of coverage across cells while fitting the count models. The models we compare use log link function, and therefore, the offsets should be log transformed too.

```r
> exposure <- log(colSums(cntmat))
> exposure
  [1] 8.815518 8.996157 9.025816 9.037771 9.062420 9.397732 8.999249 9.445412
  [9] 7.826443 9.511851 9.501367 9.132271 8.768263 9.600015 9.250042 8.235626
 [17] 9.326255 8.101678 8.439232 8.723719 9.858804 9.205830 8.708309 8.266164
 [25] 8.831128 7.811163 7.525640 8.498418 8.898775 8.716372 9.151757 9.353748
 [33] 8.915029 8.895904 8.852379 8.697012 9.550805 9.277812 9.560293 8.018625
 [41] 7.937017 8.744807 8.803875 8.259199 8.454040 7.911324 9.345046 7.434257
 [49] 8.513185 8.615771 9.757768 9.370587 7.293698 7.874359 8.672315 8.921324
 [57] 7.952967 7.960673 9.070158 9.256842 8.929170 8.144679 8.885026 9.296977
 [65] 9.161990 8.794219 9.038603 8.281724 8.206856 8.462948 8.919186 9.123038
 [73] 7.870548 7.959276 7.735433 8.036250 9.060099 9.051345 9.135293 9.806756
 [81] 7.651120 9.223257 8.799360 8.027477 9.538564 8.731336 8.835792 8.567316
 [89] 9.507032 8.770594 9.515617 7.277248 8.100465 9.211939 8.494334 8.605204
 [97] 9.034080 8.814330 8.989818 9.109857 9.642707 9.371694 8.190354 8.797700
[105] 9.128262 8.956738 8.248267 7.771067 9.630891 8.840001 8.915835 8.533460
[113] 7.687997 8.707483 9.039433 9.157783 8.098339 8.965973 9.211939 8.760610
[121] 8.556222 7.913521 8.799058 8.897135 8.092851 8.256867 7.728416 9.024854
[129] 8.022569 8.937875 8.387312 8.690474 9.481588 7.996317 7.698936 9.261794
[137] 9.538564 9.330875 9.243582 8.596928 8.944159 9.110741 7.865188 8.848509
[145] 9.139381 8.985696 8.480737 8.406708 7.813996 9.186253 8.628377 7.427739
[153] 9.163773 8.800415 7.804659 9.615739 8.962392 8.189522 8.698514 9.102310
[161] 7.980366 8.756210 7.455877 8.889997 9.748703 9.508962 8.808519 8.763897
[169] 9.059052 9.294406 7.793174 7.298445 8.411833 8.771215 8.190909 9.183894
[177] 9.115040 8.794673 7.836370 7.544332 7.252054 9.129130 7.771489 7.282074
[185] 9.444859 9.419466 9.513773 8.690138 8.046229 7.375256
```

We will pick a gene (one that we know it best fits to ZINB model) and load its UMI counts across cells into a data frame.

```r
> gg <- 4153
> gsymb[gg]
[1] "Cybb"
> y <- cntmat[gg,]
> y
  [1]  1  0  8  7  3  6  7  1  4  0  0  2  2 10  0  3  1  1  2  4  6  3  9  3  1
 [26]  4  1  2  3  8  0  0  1  6  2  6  3  0  0  1  0  4  1  3  0  2 17  1  2  5
 [51]  0  1  1  0  6  3  4  4  3  1  3  1  7  8  0  2  4  6  2  7  0  4  0  1  0
 [76]  2  0  8  2  0  1  4  9  3  0  2  2  3  0  3  0  0  4  1  4  3  2  3  2  2
[101]  8  0  5  8  2  3  0  3  0  3  6  0  4  2 11  5  1  2 12  2  4  2  6  8  6
[126]  3  3  1  0  2  1  1 10  2  2  0  2  2  1  0  1  0  3 11  0  7  2  2  2  3
[151]  7  0  1 10  4  1  4  5  1  0  1  3  3  0  5  4  9  4  2  0  4  1  2  4  1
[176]  2  1  4  2  1  0  1  0  1  4  6  8  9  4  3
```

We fit all four models. We can also add covariates to the data.frame and specify our own model formula.


```r
> gexpr <- data.frame(y, exposure)  # data.frame(y, exposure, celltype, sex)
> model_fit <- fit_count_models(gexpr, as.formula('y ~ 1 + offset(exposure)'))

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
       elpd_diff se_diff
model4    0.0       0.0
model2  -10.0       4.5
model3  -34.9       9.2
model1 -111.0      20.9

> select_model(elpd_loo, margin=2)
[1] 4
```

## How to cite

K. Choi, Y. Chen, D.A. Skelly, G.A. Churchill. “Bayesian model selection reveals biological origins of zero inflation in single-cell transcriptomics.” bioRxiv. doi: <a href="https://doi.org/10.1101/2020.03.03.974808"></a> (2020)
