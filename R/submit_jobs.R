#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param loomfile Expression quantity file (loom format)
#' @param num_chunks Number of chunks
#' @param outdir Name of the folder where results should be stored
#' @param dryrun TRUE if you only want to check the job submission commands
#' @param rfile R script to execute
#' @param scriptfile PBS job submission script
#' @param covariates Name of covariates (A col.attrs name in the input loomfile)
#' @param layer Layer name of the count to use in the loom file
#' @param nCores Number of cores to run stan fitting in parallel
#' @param seed Seed number to reproduce randomized results
#' @param gene_start Starting gene index to analyze
#' @param gene_end Ending gene index to analyze
#' @param chunk_start Starting chunk index to submit
#' @param chunk_end Ending chunk index to submit
#' @return ... None is returned
#'
submit_jobs <- function(loomfile, num_chunks, outdir, dryrun, rfile, scriptfile,
                        covariates=c(), layer=NULL, nCores=NULL, seed=NULL,
                        gene_start=NULL, gene_end=NULL, chunk_start=NULL, chunk_end=NULL) {
  if(is.null(nCores)) {
    nCores <- min(4, parallel::detectCores())
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  ds <- connect(loomfile, mode = 'r+')
  if(is.null(layer)) {
    dmat <- ds$matrix[,]
    message('[scRATE::submit_jobs] Counts from main layer will be loaded.')
  } else {
    dmat <- ds$layers[[layer]][,]
    message(sprintf('[scRATE::submit_jobs] Counts from %s layer will be loaded.', layer))
  }
  num_cells <- dim(dmat)[1]
  num_genes <- dim(dmat)[2]
  gname <- ds$row.attrs$GeneID[]
  cname <- ds$col.attrs$CellID[]
  if(length(covariates) > 0) {
    covar_list <- list()
    for (covar in covariates) {
      covar_list[[covar]] <- ds$col.attrs[[covar]][]
      message(sprintf('[scRATE::submit_jobs] %s is successfully loaded from the loom file.', covar))
    }
  } else {
    covar_list <- NULL
    message('[scRATE::submit_jobs] No covariates will be used.')
  }
  selected <- ds$row.attrs$Selected[] > 0
  ds$close_all()

  if(is.null(gene_start)) {
    gidx1 <- 1
  } else {
    gidx1 <- gene_start
  }
  if(is.null(gene_end)) {
    gidx2 <- num_genes
  } else {
    gidx2 <- gene_end
  }

  idx_gsurv <- which(selected > 0)
  idx_gsurv <- idx_gsurv[idx_gsurv >= gidx1 & idx_gsurv <= gidx2]
  num_gsurv <- length(idx_gsurv)
  message(sprintf('[scRATE::submit_jobs] %d genes (between Gene %d and %d) will be processed.', num_gsurv, gidx1, gidx2))

  if(num_chunks >= num_gsurv) {
    chunk_sz <- 1
  } else {
    chunk_sz <- num_gsurv / num_chunks
  }
  chunk_end_idx <- round(chunk_sz * 1:num_chunks)
  gene_ends <- idx_gsurv[chunk_end_idx]
  gene_starts <- gene_ends + 1
  gene_starts <- c(gidx1, gene_starts)
  gene_starts <- gene_starts[-length(gene_starts)]
  gene_ends[length(gene_ends)] <- gidx2

  dmat <- as.data.frame(t(dmat))
  rownames(dmat) <- gname
  colnames(dmat) <- cname
  csize <- as.vector(colSums(dmat))

  if (is.null(chunk_start)) {
    cidx1 <- 1
  } else {
    cidx1 <- chunk_start
  }
  if (is.null(chunk_end)) {
    cidx2 <- num_chunks
  } else {
    cidx2 <- chunk_end
    if (chunk_end > num_chunks) {
      cidx2 <- num_chunks
      message(sprintf("[scRATE::submit_jobs] There are %d chunks only, but you requested more up to %d. The last chunk index is modified accordingly.",
                  num_chunks, chunk_end))
    }
  }
  message(sprintf('[scRATE::submit_jobs] Chunk %d to %d (out of %d) will be processed.', cidx1, cidx2, num_chunks))

  for (k in cidx1:cidx2) {
    s <- gene_starts[k]
    e <- gene_ends[k]
    gsurv  <- selected[s:e]
    cntmat <- dmat[s:e,][gsurv, ]

    ifile <- file.path(outdir, sprintf('_chunk.%05d', k))
    ofile <- file.path(outdir, sprintf('_scrate.%05d', k))
    cmdstr <- sprintf('qsub -o %s -e %s -v RFILE=%s,INFILE=%s,OUTFILE=%s,CORES=%d,SEED=%d %s',
                      outdir, outdir, rfile, ifile, ofile, nCores, seed, scriptfile)
    if(!dryrun) {
      save(cntmat, csize, covar_list, file = ifile)
      message(sprintf("[scRATE::submit_jobs] Created input file: %s", ifile))
      message(sprintf("[scRATE::submit_jobs] Submitting a job: %s", cmdstr))
      system(cmdstr)
      Sys.sleep(1)
    } else {
      message(sprintf("[scRATE::submit_jobs] Will submit a job: %s", cmdstr))
    }
  }
}
