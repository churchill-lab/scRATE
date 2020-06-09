## Emacs, make this -*- mode: sh; -*-
## Adapted from rocker/r-base

FROM debian:testing

LABEL org.label-schema.license="GPLv3.0" \
      maintainer="Kwangbom (KB) Choi, Ph.D. <kb.choi@jax.org>" \
      version="0.1.1"

## Set a default user. Available via runtime flag `--user docker`
## Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
## User should also have & own a home directory (for rstudio or linked volumes to work properly).
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		fonts-texgyre \
    gcc-9-base \
    gcc \
    g++ \
    build-essential \
    libhdf5-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
  && apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## Use Debian unstable via pinning -- new style via APT::Default-Release
RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
        && echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default

ENV R_BASE_VERSION 4.0.0

RUN mkdir -p $HOME/.R/ \
    && echo "CXX14FLAGS=-O3 -march=native -mtune=native -fPIC\n" >> $HOME/.R/Makevars \
    && echo "CXX14=g++\n" >> $HOME/.R/Makevars

## Now install R and littler, and create a link for littler in /usr/local/bin
RUN apt-get update \
	&& apt-get install -t unstable -y --no-install-recommends \
		littler \
                r-cran-littler \
		r-base=${R_BASE_VERSION}-* \
		r-base-dev=${R_BASE_VERSION}-* \
		r-recommended=${R_BASE_VERSION}-* \
        r-cran-devtools \
        r-cran-rcpp \
        r-cran-rcppeigen \
        r-cran-dt \
        r-cran-bh \
        r-cran-stanheaders \
        r-cran-rcppparallel \
        r-cran-rstan \
        r-cran-rstanarm \
        r-cran-shinystan \
        r-cran-brms \
	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

RUN R --slave -e 'install.packages("rstanarm")'
RUN R --slave -e 'devtools::install_version("rstantools", version = "2.0.0", repos = "http://cran.us.r-project.org")'
RUN R --slave -e 'install.packages("brms")'
RUN R --slave -e 'devtools::install_github("mojaveazure/loomR")'
RUN R --slave -e 'devtools::install_github("churchill-lab/scRATE", dep = FALSE, build_vignettes = TRUE)'

#ENV TINI_VERSION v0.18.0
#ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
#RUN chmod +x /tini
#ENTRYPOINT ["/tini", "--"]

CMD ["/bin/bash"]
