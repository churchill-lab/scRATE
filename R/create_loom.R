#' Create a basic loom file from a tsv file that contains UMI count matrix
#'
#' @export
#' @param efile A tsv file that containes UMI count matrix (gene x cell)
#' @param outfile Name of loom file to encode single-cell gene expression
#' @param min_expressed_cells The number of non-zero count cells for a gene to be marked 'Selected' for scRATE analysis
#' @return Nothing
#'
create_loom <- function(efile, outfile, min_expressed_cells=0) {
  cntmat <- read.table(efile, row.names=1, header=TRUE, check.names=FALSE)
  message(sprintf('[scRATE::create_loom] Loaded a count matrix of %d genes x %d cells', dim(cntmat)[1], dim(cntmat)[2]))
  binmat <- cntmat > 0
  num_expr <- rowSums(binmat)
  selected <- num_expr > min_expressed_cells
  message(sprintf('[scRATE::create_loom] There were %d genes that pass the minimum expressed cell count requirement', sum(selected)))
  message(sprintf('[scRATE::create_loom] Creating a loom file'))
  loomR::create(filename = outfile,
                data = cntmat,
                feature.attrs = list(GeneID=rownames(cntmat),
                                     Selected=selected),
                cell.attrs = list(Size=colSums(cntmat)))
  message(sprintf('[scRATE::create_loom] Finished creating a loom file, %s.', outfile))
}
