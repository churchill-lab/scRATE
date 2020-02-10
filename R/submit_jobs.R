#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param loomfile Expression quantity file (loom format)
#' @param num_chunks Number of chunks
#' @param outdir Name of the folder where results should be stored
#' @param dryrun TRUE if you only want to check the job submission commands
#' @param scriptfile PBS job submission script
#' @param rfile R script to execute
#' @param layer Layer name of the count to use in the loom file
#' @param model2fit_var Variable name containing the model to fit (A row.attrs name in the input loomfile)
#' @param covariate Name of covariate (A col.attrs name in the input loomfile)
#' @param nCores Number of cores to run stan fitting in parallel
#' @param seed Seed number to reproduce randomized results
#' @param gene_start Starting gene index to analyze
#' @param gene_end Ending gene index to analyze
#' @param chunk_start Starting chunk index to submit
#' @param chunk_end Ending chunk index to submit
#' @return ... None is returned
#'
submit_jobs <- function(loomfile, num_chunks, outdir, dryrun, scriptfile, rfile,
                        layer=NULL, model2fit_var=NULL, covariate=NULL, nCores=NULL, seed=NULL,
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
    cat('[submit_jobs] Counts from main layer will be loaded.\n')
  } else {
    dmat <- ds$layers[[layer]][,]
    cat(sprintf('[submit_jobs] Counts from %s layer will be loaded.\n', layer))
  }
  num_cells <- dim(dmat)[1]
  num_genes <- dim(dmat)[2]
  gname <- ds$row.attrs$GeneID[]
  cname <- ds$col.attrs$CellID[]
  selected <- ds$row.attrs$Selected[]
  if(is.null(model2fit_var)) {
    modelID <- NULL
    cat('[submit_jobs] No specific model is set for fitting.\n')
  } else {
    modelID <- ds$row.attrs[[model2fit_var]][]
    cat(sprintf('[submit_jobs] Models specified in %s will be used for fitting.\n', model2fit_var))
  }
  if(is.null(covariate)) {
    ctype <- NULL
    cat('[submit_jobs] No covariate will be used.\n')
  } else {
    ctype <- ds$col.attrs[[covariate]][]
    cat(sprintf('[submit_jobs] %s will be used as a covariate.\n', covariate))
  }
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
  cat(sprintf('[submit_jobs] %d genes (between Gene %d and %d) will be processed.\n', num_gsurv, gidx1, gidx2))

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
      cat(sprintf("[submit_jobs] There are %d chunks only, but you requested more up to %d. The last chunk index is modified accordingly.\n",
                  num_chunks, chunk_end))
    }
  }
  cat(sprintf('[submit_jobs] Chunk %d to %d (out of %d) will be processed.\n', cidx1, cidx2, num_chunks))

  for (k in cidx1:cidx2) {
    s <- gene_starts[k]
    e <- gene_ends[k]
    cntmat <- dmat[s:e,]
    gsurv  <- selected[s:e]
    model2fit <- modelID[s:e]
    ifile <- sprintf('%s/_chunk.%05d-%05d.rds', outdir, s, e)
    ofile <- sprintf('%s/_scrate_model_fit.%05d-%05d.rds', outdir, s, e)
    cmdstr <- sprintf('qsub -o %s -e %s -v RFILE=%s,INFILE=%s,OUTFILE=%s,CORES=%d,SEED=%d %s',
                      outdir, outdir, rfile, ifile, ofile, nCores, seed, scriptfile)
    if(!dryrun) {
      save(cntmat, gsurv, csize, ctype, model2fit, file = ifile)
      cat(cmdstr, '\n')
      system(cmdstr)
      Sys.sleep(1)
    } else {
      cat(cmdstr, '\n')
    }
  }
}
