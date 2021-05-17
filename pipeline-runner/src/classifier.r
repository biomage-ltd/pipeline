################################################
## STEP 3. Classifier filter 
#################################################
#
#' @description Filters seurat object based on classifier filter using emptyDrops
#               https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_classifier)
#'          - filterSettings: slot with thresholds
#'                  - minProbabiliy: 
#'                  - filterThreshold: 
#' @export return a list with the filtered seurat object by probabilities classifier, the config and the plot values

source('utils.r')

generate_default_values_classifier <- function(seurat_obj, config) {
   
        # HARDCODE
        threshold <- 0.01
   
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
  # config$filterSettings = list(FDR=0.82, bandwidth=-1, filterThreshold=-1)
  # As of 26/3/2021
  # PlotUUID: classifierEmptyDropsPlot
  # Labels & Inputs : 
  # X-axis : cellSize, transformed with log10
  # Y-axis : FDR of emptyDrops classifier [0,1]   
  # Plot Data Schema
  # Data point :
  # { 
  #    classifierP: classifier probability,
  #    log_u: log cell size
  # }
  # Example :
  # [
  #    {"FDR": 0.994553522823595,
  #     "log_u": 4.25168687816849},
  # ...]
  print(paste("Running",task_name,"config: ",sep=" "))
  print(config)
  print(paste0("Cells per sample before filter for sample ", sample_id))
  
  print(table(seurat_obj$orig.ident, useNA="ifany"))
  #The format of the sample_id is
  # sample-WT1
  # we need to get only the last part, in order to grep the object.
  tmp_sample <- sub("sample-","",sample_id)

  import::here(map2, .from = purrr)

  # Check wheter the filter is set to true or false
  # For some reason the last children of named lists are computed as vectors, so we can't access them as recursive objects. 
  FDR <- config$filterSettings$FDR
  # Check if it is required to compute sensible values. From the function 'generate_default_values_classifier', it is expected
  # to get a list with two elements {minProbabiliy and filterThreshold}.
  if (isTRUE(config$auto)){
      FDR <- generate_default_values_classifier(seurat_obj, config)
  }
  # TODO: get flag from here: seurat_obj@tools$flag_filtered <- FALSE
  if(config$enabled) {
    # check if filter data is actually available
    if (is.null(seurat_obj@meta.data$emptyDrops_FDR)) {
      print(paste("Running",task_name,"config: ",sep=" "))
      print("Classify is enabled but has no classify data available: will dissable it: no filtering!")      
      # should this be json, i.e. "false" instead?
      config$enabled <- FALSE
      plots <- list()
    } else { # enabled and good data:
      print(paste0("Classify is enabled but classify data available: all good for filtering with FDR=", FDR))
      obj_metadata <- seurat_obj@meta.data
      # extract plotting data of original data to return to plot slot later
      barcode_names_this_sample <- rownames(obj_metadata[grep(tmp_sample, rownames(obj_metadata)),]) 
      if(length(barcode_names_this_sample)==0){
          plots <- list()
          plots[generate_plotuuid(sample_id, task_name, 0)] <- list()
          return(list(data = seurat_obj,config = config,plotData = plots)) 
      }
      sample_subset <- subset(seurat_obj, cells = barcode_names_this_sample)
      print("Info: empty-drops table of FDR threshold categories (# UMIs for a given threshold interval")
      print(table(obj_metadata$orig.ident, cut(obj_metadata$emptyDrops_FDR, breaks = c(-Inf,0,0.0001,0.01,0.1,0.5,1,Inf)), useNA="ifany"))
      
      # prevents filtering of NA FDRs if FDR=1
      ed_fdr <- sample_subset$emptyDrops_FDR
      ed_fdr[is.na(ed_fdr)] <- 1
      
      print("How many barcodes should be filtered out for this sample (#FALSE):")
      print(table(ed_fdr <= FDR))

      numis <- log10(sample_subset@meta.data$nCount_RNA)

      plot1_data <- unname(purrr::map2(ed_fdr, numis, function(x,y){c("FDR"=x,"log_u"=y)}))
      # cells: Cell names or indices
      # extract cell id that do not(!) belong to current sample (to not apply filter there)
      barcode_names_non_sample <- rownames(obj_metadata[-grep(tmp_sample, rownames(obj_metadata)),])
      # all barcodes that match threshold in the subset data
      barcode_names_keep_current_sample <-colnames(sample_subset[, ed_fdr <= FDR])
      # combine the 2:
      barcodes_to_keep <- union(barcode_names_non_sample, barcode_names_keep_current_sample)
      seurat_obj <- subset_safe(seurat_obj,barcodes_to_keep)
      # update config
      config$filterSettings$FDR <- FDR

      # Downsample plotData
      # Handle when the number of remaining cells is less than the number of cells to downsample
      num_cells_to_downsample <- downsample_plotdata(ncol(sample_subset), num_cells_to_downsample)
      print(paste('sample of size', ncol(sample_subset), 'downsampled to', num_cells_to_downsample, 'cells'))
      
      set.seed(123)
      cells_position_to_keep <- sample(1:ncol(sample_subset), num_cells_to_downsample, replace = FALSE)
      cells_position_to_keep <- sort(cells_position_to_keep)
      plot1_data <- plot1_data[cells_position_to_keep]

      #Populate plots list
      plots <-list()
      plots[generate_plotuuid(sample_id, task_name, 0)] <- list(plot1_data)
    }
  } else {
    print("filter disabled: data not filtered!")
    plots<-list()
    plots[generate_plotuuid(sample_id, task_name, 0)] <- list()
    seurat_obj <- seurat_obj
  }
  print(paste0("Cells per sample after filter for sample ", sample_id))
  print(table(seurat_obj$orig.ident, useNA="ifany"))
  # > head(plots[[1]])
  # [[1]]
  # FDR    log_u
  # 0.000000 3.554852

  # [[2]]
  # FDR        log_u
  # 0.0001120495 2.7176705030
  # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
  result <- list(
                 data = seurat_obj,
                 config = config,
                 plotData = plots
  )

  return(result)
}
