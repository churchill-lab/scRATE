% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_loom.R
\name{create_loom}
\alias{create_loom}
\title{Create a basic loom file from a tsv file that contains UMI count matrix}
\usage{
create_loom(efile, outfile, min_expressed_cells = 0)
}
\arguments{
\item{efile}{A tsv file that containes UMI count matrix (gene x cell)}

\item{outfile}{Name of loom file to encode single-cell gene expression}

\item{min_expressed_cells}{The number of non-zero count cells for a gene to be marked 'Selected' for scRATE analysis}
}
\value{
Nothing
}
\description{
Create a basic loom file from a tsv file that contains UMI count matrix
}
