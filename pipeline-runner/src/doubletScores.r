################################################
## STEP 5. Doublet score filter 
#################################################
# This is a simplest filter that looks at a threshold value for the doublet scores. 
# To separate cells with low droplet score from the ones that have a high droplet score content what makes us think that the are mistakenly considered as a single cell but they are actully two or more.  
# This can be a useful first guess. The settings for such a filter can also contain a simple "probabilityThreshold" setting. 

source('utils.r')

# The most uses values in doublet scores reporting in the scrublet paper [1] are around 0.25. There are not too much literature about how to compute
# a threshold. For now, we will offer two methods:
# --> Absolute threshold: In order to be not too extrictive the threshold is set to 0.25
generate_default_values_doubletScores <- function(seurat_obj, config) {

    # HARDCODE
    threshold <- 0.25

    return(threshold)
}

#' @description Filters seurat object based on doubletScores scores
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_doubletScores)
#'          - filterSettings: slot with thresholds
#'                  - probabilityThreshold: Float. cut-off for the maximun probability scores that could have a cell. The doubletScores scores have been computed
#'                  through "scrublet" [1]
#'                  - binStep: Float. Bin size for the histogram
#' @export return a list with the filtered seurat object by doublet score, the config and the plot values

task <- function(seurat_obj, config, task_name, sample_id, num_cells_to_downsample = 6000){
    print(paste("Running",task_name,sep=" "))
    print("Config:")
    print(config)
    print(paste0("Cells per sample before filter for sample ", sample_id))
    print(table(seurat_obj$orig.ident))
    #The format of the sample_id is
    # sample-WT1
    # we need to get only the last part, in order to grep the object.
    tmp_sample <- sub("sample-","",sample_id)

    # Check if the experiment has doubletScores
    if (!"doublet_scores"%in%colnames(seurat_obj@meta.data)){
        message("Warning! No doubletScores scores has been computed for this experiment!")
        return(seurat_obj)
    }
    probabilityThreshold <- config$filterSettings[["probabilityThreshold"]]

    # Check if it is required to compute sensible values. From the function 'generate_default_values_doubletScores', it is expected
    # to get a value --> probabilityThreshold.
    if (exists('auto', where=config)){
        if (as.logical(toupper(config$auto)))
            probabilityThreshold <- generate_default_values_doubletScores(seurat_obj, config)
    }

    # extract plotting data of original data to return to plot slot later
    obj_metadata <- seurat_obj@meta.data
    barcode_names_this_sample <- rownames(obj_metadata[grep(tmp_sample, rownames(obj_metadata)),]) 
    if(length(barcode_names_this_sample)==0){
        plots <- list()
        plots[generate_plotuuid(sample_id, task_name, 0)] <- list()
        return(list(data = seurat_obj,config = config,plotData = plots)) 
    }
    sample_subset <- subset(seurat_obj, cells = barcode_names_this_sample)

    # Check wheter the filter is set to true or false
    if (as.logical(toupper(config$enabled))){
        # extract cell id that do not(!) belong to current sample (to not apply filter there)
        barcode_names_non_sample <- rownames(obj_metadata[-grep(tmp_sample, rownames(obj_metadata)),]) 
        # all barcodes that match threshold in the subset data
        barcode_names_keep_current_sample <-rownames(sample_subset@meta.data[sample_subset$doublet_scores <= probabilityThreshold,])
        # combine the 2:
        barcodes_to_keep <- union(barcode_names_non_sample, barcode_names_keep_current_sample)
        # Information regarding doublet score is pre-computed during the 'data-ingest'. 
        seurat_obj.filtered <- subset_safe(seurat_obj,barcodes_to_keep)
    }
    else{
        seurat_obj.filtered <- seurat_obj
    }
        
    # update config
    config$filterSettings$probabilityThreshold <- probabilityThreshold

    plot1_data <- lapply(unname(sample_subset$doublet_scores),function(x) {c("doubletP"=x)})


    # Downsample plotData
    # Handle when the number of remaining cells is less than the number of cells to downsample
    num_cells_to_downsample <- downsample_plotdata(ncol(sample_subset), num_cells_to_downsample)
    print(paste('sample of size', ncol(sample_subset), 'downsampled to', num_cells_to_downsample, 'cells'))
    
    set.seed(123)
    cells_position_to_keep <- sample(1:ncol(sample_subset), num_cells_to_downsample, replace = FALSE)
    cells_position_to_keep <- sort(cells_position_to_keep)
    plot1_data <- plot1_data[cells_position_to_keep]

    plots <-list()

    # plot 1: histgram of doublet scores
    #              [0.161,              0.198,              0.284,  ...]
    plots[generate_plotuuid(sample_id, task_name, 0)] <- list(plot1_data)
    
    print(paste0("Cells per sample after filter for sample ", sample_id))
    print(table(seurat_obj.filtered$orig.ident))

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = seurat_obj.filtered,
        config = config,
        plotData = plots
    )
    return(result)

}


# [1] Wolock SL, Lopez R, Klein AM. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. 
# Cell Syst. 2019 Apr 24;8(4):281-291.e9. doi: 10.1016/j.cels.2018.11.005. Epub 2019 Apr 3. PMID: 30954476; PMCID: PMC6625319.
