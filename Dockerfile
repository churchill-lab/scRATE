## Emacs, make this -*- mode: sh; -*-
## Adapted from rocker/r-base

FROM r-base

LABEL org.label-schema.license="GPLv3.0" \
      maintainer="Kwangbom (KB) Choi, Ph.D. <kb.choi@jax.org>" \
      version="0.1.2"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libhdf5-dev \
        libv8-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p $HOME/.R/ \
    && echo "CXX14FLAGS=-O3 -march=native -mtune=native -fPIC\n" >> $HOME/.R/Makevars

RUN R --slave -e 'install.packages("Rcpp")'
RUN R --slave -e 'install.packages("RcppParallel")'
RUN R --slave -e 'install.packages("rstan")'
RUN R --slave -e 'install.packages("rstanarm")'
RUN R --slave -e 'install.packages("rstantools")'
RUN R --slave -e 'install.packages("brms")'
RUN R --slave -e 'install.packages("devtools")'
RUN R --slave -e 'devtools::install_github("mojaveazure/loomR")'
RUN R --slave -e 'devtools::install_github("churchill-lab/scRATE", dep = FALSE, build_vignettes = TRUE)'

CMD ["/bin/bash"]
