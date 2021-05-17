################################################
## STEP 2. Mitochondrial content filter 
#################################################
# This is a simplest filter that looks at a threshold value for the mitochondrial content.
# To separate cells with low MT-content from the ones that have a high MT-content what makes us think that are dead.  
# This can be a useful first guess. The settings for such a filter can also contain a simple "probabilityThreshold" setting. 

source('utils.r')

# The most uses values in MT-content are between [0.1, 0.2]. There are not too much literature about how to compute
# a threshold. For now, we will offer two methods:
# --> Absolute threshold: In order to be not too extrictive the threshold is set to 0.1
generate_default_values_mitochondrialContent <- function(seurat_obj, config) {
   
   if (config$filterSettings$method == "absolute_threshold")
        # HARDCODE
        threshold <- 0.1
   
    return(threshold)
}

#' @description Filters seurat object based on mitochondrialContent
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_mitochondrialContent)
#'          - filterSettings: slot with thresholds
#'                  - method: String. Method to be used {absolute_threshold} 
#'                  - methodSettings: List with the method as key and contain all the filterSettings for this specific method. 
#'                          * absolute_threshold: based on a cut-off threshold
#'                                  - maxFraction: Float. maximun pct MT-content that we considere for a alive cell
#'                                  - binStep: Float. Bin size for the histogram
#'                          * we are supposed to add more methods ....
#' @export return a list with the filtered seurat object by mitochondrial content, the config and the plot values

task <- function(seurat_obj, config, task_name, sample_id, num_cells_to_downsample = 6000){
    print(paste("Running",task_name,sep=" "))
    print("Config:")
    print(config)
    print(paste0("Cells per sample before filter for sample ", sample_id))
    print(table(seurat_obj$orig.ident))
    tmp_sample <- sub("sample-","",sample_id)
    # Check if the experiment has MT-content
    if (!"percent.mt"%in%colnames(seurat_obj@meta.data)){
        #This return value breaks the step, no config data is returned!!!!
        message("Warning! No MT-content has been computed for this experiment!")
        return(seurat_obj)
    }
    #The absolute threshold object is an atomic vector, not a list!
    maxFraction <- config$filterSettings$methodSettings[[config$filterSettings$method]]$maxFraction
    # Check if it is required to compute sensible values. From the function 'generate_default_values_doubletScores', it is expected
    # to get a list with only one element --> maxFraction.
    #Computing the subset for the current sample
    obj_metadata <- seurat_obj@meta.data
    barcode_names_this_sample <- rownames(obj_metadata[grep(tmp_sample, rownames(obj_metadata)),]) 
    if (length(barcode_names_this_sample)==0){
        plots <- list()
        plots[generate_plotuuid(sample_id, task_name, 0)] <- list()
        plots[generate_plotuuid(sample_id, task_name, 1)] <- list()
        return(list(data = seurat_obj,config = config,plotData = plots)) 
    }
    sample_subset <- subset(seurat_obj, cells = barcode_names_this_sample)

    if (exists('auto', where=config)){
        if (as.logical(toupper(config$auto)))
        maxFraction <- generate_default_values_mitochondrialContent(sample_subset, config)
    }
    # Check whether the filter is set to true or false
    if (as.logical(toupper(config$enabled))){
        # extract cell id that do not(!) belong to current sample (to not apply filter there)
        barcode_names_non_sample <- rownames(obj_metadata[-grep(tmp_sample, rownames(obj_metadata)),]) 
        # all barcodes that match threshold in the subset
        barcode_names_keep_current_sample <-rownames(sample_subset@meta.data[sample_subset@meta.data$percent.mt <= maxFraction*100,])
        # combine the 2:
        barcodes_to_keep <- union(barcode_names_non_sample, barcode_names_keep_current_sample)
        # TODO: try(): make sure  barcodes_to_keep is not empty
        # Information regarding MT-content is pre-computed during the 'data-ingest'. 
        seurat_obj.filtered <- subset_safe(seurat_obj, barcodes_to_keep)
    } else {
        seurat_obj.filtered <- seurat_obj
    }
    plot1_data <- lapply(unname(sample_subset$percent.mt),function(x) {c("fracMito"=x)})
    plot2_data <- unname(purrr::map2(sample_subset$percent.mt,sample_subset$nCount_RNA,function(x,y){c("cellSize"=y,"fracMito"=x)}))
    # update config
    config$filterSettings$methodSettings[[config$filterSettings$method]]$maxFraction <- maxFraction
    
    # Downsample plotData
    # Handle when the number of remaining cells is less than the number of cells to downsample
    num_cells_to_downsample <- downsample_plotdata(ncol(sample_subset), num_cells_to_downsample)
    print(paste('sample of size', ncol(sample_subset), 'downsampled to', num_cells_to_downsample, 'cells'))
      
    set.seed(123)
    cells_position_to_keep <- sample(1:ncol(sample_subset), num_cells_to_downsample, replace = FALSE)
    cells_position_to_keep <- sort(cells_position_to_keep)
    plot1_data <- plot1_data[cells_position_to_keep]
    plot2_data <- plot2_data[cells_position_to_keep]

    plots <- list()

    # plot 1: histgram of MT-content
    # AAACCCAAGCGCCCAT-1 AAACCCAAGGTTCCGC-1 AAACCCACAGAGTTGG-1
    #              0.161              0.198              0.284  ...
    plots[generate_plotuuid(sample_id, task_name, 0)] <- list(plot1_data)

    # plot 2: There are two alternavitive:
    #           - Scatter plot with UMIs in the x-axis and MT-content in the y-axis
    #           --> code: plot2 = list(u=seurat_obj$nCount_RNA.mt, "MT-content" = seurat_obj$percent.mt)
    #           - Barplot representing in the x-axis the log10(UMIs) and in the y-axis the MT-content. This option is the one 
    #           that is shown in the mockup.
    #           --> code: plot2 = list(log_10_UMIs=log10(seurat_obj$nCount_RNA), MT_content =mseurat_obj$percent.mt)
    # We have decided to use the scatter plot, but I temporaly leave the other option in the comments. 
    # Q: Should we return from the R side the cells that are going to be removed? For this plot it is interesting to color the
    # cells that are going to be excluded. 
    plots[generate_plotuuid(sample_id, task_name, 1)] <- list(plot2_data)

    # some tests:
    print(paste0("Cells per sample after filter for sample ", sample_id))
    print(table(seurat_obj.filtered$orig.ident))

    result <- list( 
        data = seurat_obj.filtered, # seurat_obj filter
        config = config,
        plotData = plots
    )
    return(result)
}
