#' Divide the entire data into many chunks to "run_model_comparison" on a cluster in parallel
#'
#' @export
#' @param loomfile Expression quantity file (loom format).
#' @param num_chunks Number of chunks.
#' @param outdir Name of the folder where results should be stored.
#' @param dryrun TRUE if you only want to check the job submission commands.
#' @param covariate Name of covariate (A col.attrs name in the input loomfile.)
#' @param layer A layer name of the count to use in the loom file.
#' @param nCores Number of cores to run stan fitting in parallel
#' @param seed Seed number to reproduce randomized results
#' @param gene_start Starting gene index to analyze.
#' @param gene_end Ending gene index to analyze.
#' @param chunk_start Starting chunk index to submit.
#' @param chunk_end Ending chunk index to submit.
#' @return ... None is returned.
#'
prepare_job_array <- function(loomfile, num_chunks, outdir, dryrun,
                              covariate=c(), layer=NULL, nCores=NULL, seed=NULL,
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
  covar_list <- list()
  if(length(covariate) > 0) {
    for (covar in covariate) {
      covar_list[[cover]] <- ds$col.attrs[[covar]][]
      cat(sprintf('[prepare_job_array] %s is successfully loaded from the loom file.\n', covar))
    }
  } else {
    cat('[prepare_job_array] No covariate will be used.\n')
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
      gsurv  <- selected[s:e]
      # cntmat <- dmat[s:e,]
      cntmat <- dmat[s:e,][gsurv, ]
      outfile <- file.path(outdir, sprintf('_chunk.%05d', k))
      save(cntmat, csize, covar_list, file = outfile)
      cat(sprintf("[prepare_job_array] Created input file: %s\n", outfile))
    }
    sh_file_slurm <- file.path(outdir, 'run_subjobs.sh.slurm')
    cat('#!/bin/bash\n', file=sh_file_slurm)
    cat('#SBATCH --qos=batch\n', file=sh_file_slurm, append=TRUE)
    cat('#SBATCH --partition=compute\n', file=sh_file_slurm, append=TRUE)
    cat('#SBATCH --job-name=scRATE\n', file=sh_file_slurm, append=TRUE)
    cat('#SBATCH --nodes=1\n', file=sh_file_slurm, append=TRUE)
    cat('#SBATCH --ntasks=1\n', file=sh_file_slurm, append=TRUE)
    cat('#SBATCH --cpus-per-task=4\n', file=sh_file_slurm, append=TRUE)
    cat('#SBATCH --mem=16gb\n', file=sh_file_slurm, append=TRUE)
    cat('#SBATCH --time=23:59:59\n', file=sh_file_slurm, append=TRUE)
    cat(sprintf('#SBATCH --array=1-%d\n\n', num_chunks), file=sh_file_slurm, append=TRUE)
    cat('module load singularity\n\n', file=sh_file_slurm, append=TRUE)
    cat('# Run your R script that uses scRATE singularity container\n', file=sh_file_slurm, append=TRUE)
    cat('ARRAY_ID=`printf %05d $SLURM_ARRAY_TASK_ID`\n', file=sh_file_slurm, append=TRUE)
    cat('singularity run --app Rscript ${CONTAINER} ${RFILE} _chunk.${ARRAY_ID} _scrate.${ARRAY_ID} ${CORES} ${SEED}\n', file=sh_file_slurm, append=TRUE)
    Sys.chmod(sh_file_slurm, '0755')
    # Script for PBS/Torque
    sh_file_pbs <- file.path(outdir, 'run_subjobs.sh.pbs')
    cat('#!/bin/bash\n', file=sh_file_pbs)
    cat('#PBS -l nodes=1:ppn=4\n', file=sh_file_pbs, append=TRUE)
    cat('#PBS -l mem=16gb\n', file=sh_file_pbs, append=TRUE)
    cat('#PBS -l walltime=23:59:59\n', file=sh_file_pbs, append=TRUE)
    cat(sprintf('#PBS -t 1-%d\n\n', num_chunks), file=sh_file_pbs, append=TRUE)
    cat('cd $PBS_O_WORKDIR\n\n', file=sh_file_pbs, append=TRUE)
    cat('# Add environment modules for your R here, for example ...\n', file=sh_file_pbs, append=TRUE)
    cat('module load hdf5/1.8.14\n', file=sh_file_pbs, append=TRUE)
    cat('module load R/3.5.1\n\n', file=sh_file_pbs, append=TRUE)
    cat('# Run your R script that uses scRATE library\n', file=sh_file_pbs, append=TRUE)
    cat('ARRAY_ID=`printf %05d $PBS_ARRAYID`\n', file=sh_file_pbs, append=TRUE)
    cat('Rscript ${RFILE} _chunk.${ARRAY_ID} _scrate.${ARRAY_ID} ${CORES} ${SEED}\n', file=sh_file_pbs, append=TRUE)
    Sys.chmod(sh_file_pbs, '0755')
    cat(sprintf("[prepare_job_array] Generated bash script, %s, for submitting array jobs. Modify the file if needed.\n", sh_file_pbs))
  } else {
    for (k in cidx1:cidx2) {
      outfile <- file.path(outdir, sprintf('_chunk.%05d', k))
      cat(sprintf("[prepare_job_array::dryrun] Will created input file: %s\n", outfile))
    }
  }
}
