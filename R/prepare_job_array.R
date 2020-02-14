#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param loomfile Expression quantity file (loom format).
#' @param num_chunks Number of chunks.
#' @param outdir Name of the folder where results should be stored.
#' @param dryrun TRUE if you only want to check the job submission commands.
#' @param layer Layer name of the count to use in the loom file.
#' @param covariate Name of covariate (A col.attrs name in the input loomfile.)
#' @param nCores Number of cores to run stan fitting in parallel
#' @param seed Seed number to reproduce randomized results
#' @param hpc Cluster HPC system 'slurm' (default) or 'pbs'
#' @param gene_start Starting gene index to analyze.
#' @param gene_end Ending gene index to analyze.
#' @param chunk_start Starting chunk index to submit.
#' @param chunk_end Ending chunk index to submit.
#' @return ... None is returned.
#'
prepare_job_array <- function(loomfile, num_chunks, outdir, dryrun,
                              layer=NULL, covariate=NULL, nCores=NULL, seed=NULL, hpc=c('slurm', 'pbs'),
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
    cat('[prepare_job_array] Counts from main layer will be loaded.\n')
  } else {
    dmat <- ds$layers[[layer]][,]
    cat(sprintf('[prepare_job_array] Counts from %s layer will be loaded.\n', layer))
  }
  num_cells <- dim(dmat)[1]
  num_genes <- dim(dmat)[2]
  gname <- ds$row.attrs$GeneID[]
  cname <- ds$col.attrs$CellID[]
  if(is.null(covariate)) {
    ctype <- NULL
    cat('[prepare_job_array] No covariate will be used.\n')
  } else {
    ctype <- ds$col.attrs[[covariate]][]
    cat(sprintf('[prepare_job_array] %s will be used as a covariate.\n', covariate))
  }
  selected <- ds$row.attrs$Selected[]
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
  cat(sprintf('[prepare_job_array] %d genes (between Gene %d and %d) will be processed.\n', num_gsurv, gidx1, gidx2))

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
      cat(sprintf("[prepare_job_array] There are %d chunks only, but you requested more up to %d. The last chunk index is modified accordingly.\n",
                  num_chunks, chunk_end))
    }
  }
  cat(sprintf('[prepare_job_array] Chunk %d to %d (out of %d) will be created.\n', cidx1, cidx2, num_chunks))

  if(!dryrun) {
    for (k in cidx1:cidx2) {
      s <- gene_starts[k]
      e <- gene_ends[k]
      cntmat <- dmat[s:e,]
      gsurv  <- selected[s:e]
      outfile <- file.path(outdir, sprintf('_chunk.%05d', k))
      save(cntmat, gsurv, csize, ctype, file = outfile)
      cat(sprintf("[prepare_job_array] Created input file: %s\n", outfile))
    }
    sh_file <- file.path(outdir, 'run_subjobs.sh')
    if (hpc == 'slurm') {
      cat('#!/bin/bash\n', file=sh_file)
      cat('#SBATCH --qos=batch', file=sh_file, append=TRUE)
      cat('#SBATCH --partition=compute', file=sh_file, append=TRUE)
      cat('#SBATCH --job-name=scRATE', file=sh_file, append=TRUE)
      cat('#SBATCH --nodes=1', file=sh_file, append=TRUE)
      cat('#SBATCH --ntasks=1', file=sh_file, append=TRUE)
      cat('#SBATCH --cpus-per-task=4', file=sh_file, append=TRUE)
      cat('#SBATCH --mem=16gb', file=sh_file, append=TRUE)
      cat('#SBATCH --time=23:59:59', file=sh_file, append=TRUE)
      cat(sprintf('#SBATCH --array=1-%d', num_chunks), file=sh_file, append=TRUE)
      cat('module load singularity', file=sh_file, append=TRUE)
      cat('ARRAY_ID=`printf %05d $SLURM_ARRAY_TASK_ID`', file=sh_file, append=TRUE)
      cat('singularity run --app Rscript ${CONTAINER} ${RFILE} ${INFILE}.${ARRAY_ID} ${OUTFILE}.${ARRAY_ID} ${CORES} ${SEED}', file=sh_file, append=TRUE)
    } else if (hpc == 'pbs') {
      cat('#!/bin/bash\n', file=sh_file)
      cat('#PBS -l nodes=1:ppn=4,mem=16gb,walltime=23:59:59\n\n', file=sh_file, append=TRUE)
      cat('cd $PBS_O_WORKDIR\n\n', file=sh_file, append=TRUE)
      cat('# Add environment modules here, for example ...\n', file=sh_file, append=TRUE)
      cat('module load hdf5/1.8.14\n', file=sh_file, append=TRUE)
      cat('module load R/3.5.1\n\n', file=sh_file, append=TRUE)
      cat('# Run R script\n', file=sh_file, append=TRUE)
      cat('ARRAY_ID=`printf %05d $PBS_ARRAYID`\n', file=sh_file, append=TRUE)
      cat('Rscript ${RFILE} _chunk.${ARRAY_ID} _scrate_elpd_loo.${ARRAY_ID} ${CORES} ${SEED}\n', file=sh_file, append=TRUE)
    }
    Sys.chmod(sh_file, '0755')
    cat(sprintf("[prepare_job_array] Generated bash script, %s, for submitting array jobs. Modify the file if needed.\n", sh_file))
  } else {
    for (k in cidx1:cidx2) {
      outfile <- file.path(outdir, sprintf('_chunk.%05d', k))
      cat(sprintf("[prepare_job_array::dryrun] Will created input file: %s\n", outfile))
    }
  }
}
