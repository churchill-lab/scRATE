#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param fit_list A list of model fitting results
#' @param margin A multiplier for standard deviation (SD) in leave-one-out ELPD for calling models
#' @param loomfile A expression quantity file (loom format)
#' @param attr_name Name of the row attribute in loomfile for storing best model calls
#' @param verbose Whether to print out overall model calls
#' @return modelcall Models called for genes. Stored in the input loomfile if provided.
#'
perform_model_selection <- function(fit_list, margin=2, loomfile=NULL, attr_name=NULL, verbose=FALSE) {
                                     
  gsurv <- names(fit_list)
  modelcall <- c()
  for (g in gsurv) {
    modelcall <- c(modelcall, select_model(fit_list[[g]][['elpd_loo']], margin=margin))
  }

  if(verbose==TRUE) {
    cat(sprintf('%5d Poisson genes\n', sum(modelcall==1)))
    cat(sprintf('%5d Negative Binomial genes\n', sum(modelcall==2)))
    cat(sprintf('%5d Zero-Inflated Poisson genes\n', sum(modelcall==3)))
    cat(sprintf('%5d Zero-Inflated Neg. Binomial genes\n', sum(modelcall==4)))
    cat(sprintf('%5d Total\n', length(fit_list)))
  }

  if(is.null(loomfile)) {
    return(modelcall)
  } else {
    ds <- connect(loomfile, mode='r+')
    gname <- ds$row.attrs$GeneID[]
    gene_idx <- match(gsurv, gname)
    modelcall_fullsize <- rep(0, length(gname))
    modelcall_fullsize[gene_idx] <- modelcall
    ra <- vector(mode="list", length=1)
    if(is.null(attr_name)) {
      attr_name <- sprintf('ModelCall_SD%.1f', margin)
    }
    names(ra) <- attr_name
    ra[[attr_name]] <- modelcall_fullsize
    ds$add.row.attribute(ra, overwrite=TRUE)
    ds$close_all()
    return(modelcall_fullsize)
  }

}
