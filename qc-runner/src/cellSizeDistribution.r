
################################################
## STEP 1. Cell size distribution filter
#################################################
# cell_size_distribution_filter function 
# This is a simplest filter that looks at the shape of the cell size (# of UMIs per cell) distribution and looks for some local minima, minimum of second derivative, etc. 
# To separate the large cell barcodes that correspond to real cells from the tail containing 'empty droplets'. 
# This can be a useful first guess. The settings for such a filter can also contain a simple "min cell size" setting. 
#
#' @description Filters seurat object based on cell size distribution
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_cellSizeDistribution)
#'          - filterSettings: slot with thresholds
#'                  - minCellSize: Integer. Threshold that contain the minimun number of UMIs per cell
#'                  - binStep: Integer. Bin size for the histogram
#' @export return a list with the filtered seurat object by cell size ditribution, the config and the plot values

source('utils.r')

# CalculateBarcodeInflections calculates an adaptive inflection point ("knee")
# of the barcode distribution for each sample group. This is
# useful for determining a threshold for removing low-quality
# samples.
generate_default_values_cellSizeDistribution <- function(seurat_obj, config) {
  import::here(CalculateBarcodeInflections, .from = Seurat)
  import::here(Tool, .from = Seurat)
  # Q: should we check for precalculated values? e.g.:
  # is.null(Tool(seurat_obj, slot = "CalculateBarcodeInflections"))
  # Returns Seurat object with a new list in the `tools` slot,
  # `CalculateBarcodeInflections` including inflection point calculatio
  seurat_obj_tmp <- CalculateBarcodeInflections(
                                        object = seurat_obj,
                                        barcode.column = "nCount_RNA",
                                        group.column = "orig.ident",
                                        # [HARDCODED]
                                        threshold.low = 1e2,
                                        threshold.high = NULL
  )
  # extracting the inflection point value which can serve as minCellSize threshold
  # all other side effects to the scdate object will be discarded
  # [TODO unittest]: this calculation assumes a single sample, i.e. only one group in orig.ident!
  # otherwise it will return multiple values for each group!
  # returned is both the rank(s) as well as inflection point
  # orig.ident          nCount_RNA  rank
  # SeuratProject       1106        10722
  sample_subset <- Tool(seurat_obj_tmp, slot = "CalculateBarcodeInflections")$inflection_points
  # extracting only inflection point(s)
  return(sample_subset$nCount_RNA)
}

task <- function(seurat_obj, config, task_name, sample_id) {
    print(paste("Running",task_name,sep=" "))
    print("Config:")
    print(config)
    print(paste0("Cells per sample before filter for sample ", sample_id))
    print(table(seurat_obj$orig.ident))
    #The format of the sample_id is
    # sample-WT1
    # we need to get only the last part, in order to grep the object.
    tmp_sample <- sub("sample-","",sample_id)

    import::here(map2, .from = purrr)
    minCellSize <- as.numeric(config$filterSettings$minCellSize)

    # extract plotting data of original data to return to plot slot later
    obj_metadata <- seurat_obj@meta.data
    barcode_names_this_sample <- rownames(obj_metadata[grep(tmp_sample, rownames(obj_metadata)),]) 
    
    if(length(barcode_names_this_sample)==0){
        plots <- list()
        plots[generate_plotuuid(sample_id, task_name, 0)] <- list()
        plots[generate_plotuuid(sample_id, task_name, 1)] <- list()
        return(list(data = seurat_obj,config = config,plotData = plots)) 
    }
    sample_subset <- subset(seurat_obj, cells = barcode_names_this_sample)
    plot1_data  <- sort(sample_subset$nCount_RNA, decreasing = TRUE)
    # plot 1: histgram of UMIs, hence input is all UMI values, e.g.
    # AAACCCAAGCGCCCAT-1 AAACCCAAGGTTCCGC-1 AAACCCACAGAGTTGG-1
    #               2204              20090               5884  ..   
    #    u    u    u    u    u    u
    # 3483 6019 3892 3729 4734 3244
    
    plot2_data <- seq_along(plot1_data)
    plot2_data <- unname(map2(plot1_data,plot2_data,function(x,y){c("u"=x,"rank"=y)}))
    plot1_data <- lapply(unname(plot1_data),function(x) {c("u"=x)})

    # Check if it is required to compute sensible values. From the function 'generate_default_values_cellSizeDistribution', it is expected
    # to get a list with two elements {minCellSize and binStep}   
    if (exists('auto', where=config)){
        if (as.logical(toupper(config$auto)))
            # config not really needed for this one (maybe later for threshold.low/high):
            minCellSize <- generate_default_values_cellSizeDistribution(sample_subset, config)
    }
    if (as.logical(toupper(config$enabled))) {
        # Information regarding number of UMIs per cells is pre-computed during the 'CreateSeuratObject' function. 
        # this used to be the filter below which applies for all samples
        # we had to change it for the current version where we also have to 
        # keep all barcodes for all samples (filter doesn't apply there)
        # once we ensure the input is only one sample, we can revert to the line below:
        # seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= minCellSize)
        # subset(x, subset, cells = NULL, features = NULL, idents = NULL, ...)
        # cells: Cell names or indices
        # extract cell id that do not(!) belong to current sample (to not apply filter there)
        barcode_names_non_sample <- rownames(obj_metadata[-grep(tmp_sample, rownames(obj_metadata)),]) 
        # all barcodes that match threshold in the subset data
        barcode_names_keep_current_sample <-rownames(sample_subset@meta.data[sample_subset@meta.data$nCount_RNA >= minCellSize,])
        # combine the 2:
        barcodes_to_keep <- union(barcode_names_non_sample, barcode_names_keep_current_sample)
        
        seurat_obj.filtered <- subset_safe(seurat_obj,barcodes_to_keep)
    }else{
        seurat_obj.filtered <- seurat_obj
    }

    print(paste0("Cells per sample after filter for sample ", sample_id))
    print(table(seurat_obj.filtered$orig.ident))

    # update config
    config$filterSettings$minCellSize <- minCellSize

    #Populate plots list
    plots <-list()
    plots[generate_plotuuid(sample_id, task_name, 0)] <- list(plot1_data)
    plots[generate_plotuuid(sample_id, task_name, 1)] <- list(plot2_data)

    result <- list(
        data = seurat_obj.filtered,
        config = config,
        plotData = plots
    )

    return(result)
}
