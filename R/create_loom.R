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
  binmat <- cntmat > 0
  num_expr <- rowSums(binmat)
  loomR::create(filename = outfile,
                data = cntmat,
                feature.attrs = list(GeneID=rownames(cntmat),
                                     Selected=num_expr > min_expressed_cells),
                cell.attrs = list(Size=colSums(cntmat)))
}
