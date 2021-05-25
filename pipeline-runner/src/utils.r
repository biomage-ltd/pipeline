#
# Generate GUI plots and data uuid as required by the UI
# 
generate_gui_uuid <- function(sample_uuid, task_name, item_idx) {

  if(sample_uuid != "") {
    return(paste(sample_uuid, task_name, item_idx, sep="-"))
  }

  return(paste(task_name, item_idx, sep="-"))

}

#
# Subset safe allows us to attempt to subset a seurat object even with an empty list
# this function exists for the case when the user over-filters the object, and we need to return something
# that'd allow the user to realize that they are filtering all the cells, while maintaining certain seurat functionality.
# it's a questionable function and it should be questioned.
#
# IN seurat_obj: object to filter
# IN cells: cell barcodes to subset the object with
# 
subset_safe <- function(seurat_obj,cells){
  if(length(cells)>0){
    return(subset(seurat_obj, cells = cells))
  }else{
    return(subset(seurat_obj,cells = colnames(seurat_obj)[1]))
  }
}


#
# down sample plots
# 

downsample_plotdata <- function(ncol_sample, max_number_of_cells) {
  return(min(max_number_of_cells, ncol_sample))
}


handle_debug <- function(scdata, config, task_name, sample_id, debug_config) {
    is_debug <- debug_config$step %in% c(task_name, 'all')
    
    if (is_debug) {
        # variable names used by functions
        seurat_obj <- scdata
        num_cells_to_downsample <- 6000    
        sample_str <- ifelse(sample_id == '', '', paste0('_', sample_id))
        fname <- paste0(task_name, sample_str, '.RData')
        fpath_cont <- file.path('/debug', fname)
        fpath_host <- file.path(debug_config$path, fname)
        message(sprintf('⚠️ DEBUG_STEP = %s. Saving arguments.', task_name))
        save(seurat_obj, config, task_name, sample_id, num_cells_to_downsample, file = fpath_cont)
        message(sprintf("⚠️ RUN load('%s') to restore environment.", fpath_host))
    }
}


#' Calculates statistics before/after filter step
#'
#' @param seurat_obj \code{SeuratObject}
#' @param tmp_sample sample name in \code{seurat_obj@orig.ident} to compute statistics for
#'
#' @return list with \itemize{
#'   \item{"num_cells"}{Number of cells in sample}
#'   \item{"total_genes"}{Number of detected genes in sample}
#'   \item{"median_genes"}{Median number of genes detected per cell}
#'   \item{"median_umis"}{Median number of counts per cell}
#' }
#' 
calc_filter_stats <- function(seurat_obj, tmp_sample) {
  
  # subset to current sample
  scdata <- seurat_obj[, seurat_obj$orig.ident == tmp_sample]
  
  # number of counts per gene
  ncount <- Matrix::rowSums(scdata[['RNA']]@counts)
  
  list(
    num_cells = ncol(scdata),
    total_genes = sum(ncount > 0),
    median_genes = median(scdata$nFeature_RNA),
    median_umis = median(scdata$nCount_RNA)
  )
}
