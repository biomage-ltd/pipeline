################################################
## 2-1_Compute-metrics_emptyDrops.r
##  - Compute emptyDrops per sample
##  - Save an rds with emptyDrops result per sample
################################################


suppressWarnings(library(Seurat))
suppressWarnings(library(DropletUtils))

# compute_emptydrops function
#' @description Save the result of emptyDrops per sample. 
#' @param scdata Raw sparse matrix with the counts for one sample.
#' @param sample_name Name of the sample that we are preparing.
#'
#' @return 
compute_emptydrops <- function(scdata, sample_name) {
  
  message("Sample --> ", sample_name, "...")

  # [HARDCODE]
  threshold_emptydrops <- 100
  # threshold_emptydrops <- 2000 # use this to simulate unfiltered data
  # automatically checking if data contains enough (50) barcodes that 
  # are confidently considered to be empty to serve as training set
  # not necessary as of now: sce <- as.SingleCellExperiment(scdata)
  nr_barcodes_ambient <- sum(colSums(scdata) < threshold_emptydrops) 
  range(colSums(scdata) ) 
  if (nr_barcodes_ambient < 50 ) {
    flag_filtered <- TRUE
    message("Found not enough [", nr_barcodes_ambient, "] ambient barcodes,
          so data is considered to be pre-filtered")
    message("Smallest nr of UMI for any barcode: ", min(colSums(scdata)))
    message("EmptyDrops filter will not be applied")
  } else {
    flag_filtered <- FALSE
    # if non-filtered, run emptyDrops and save results in file
    message("Found enough [", nr_barcodes_ambient, "] ambient barcodes, 
          so data is considered to be not filtered")
    emptydrops_out <- emptyDrops(scdata, lower = threshold_emptydrops) 

    file_path <- paste("/output/pre-emptydrops-", sample_name, ".rds", sep = "")
    saveRDS(emptydrops_out, file = file_path, compress = FALSE)
  }
}

task <- function(input,pipeline_config){
  message("reloading old matrices...")
  scdata_list <- readRDS("/output/pre-doublet-scdata_list.rds")

  message("Loading configuration...")
  config <- RJSONIO::fromJSON("/input/meta.json")

  # Check which samples have been selected. Otherwiser we are going to use all of them. 
  if (length(config$samples)>0){
      samples <- config$samples
  }else{
      samples <- names(scdata_list)
  }

  scdata_list <- scdata_list[samples]

  message("calculating probability of barcodes being background noise...")
  for (sample_name in names(scdata_list)) {
    compute_emptydrops(scdata_list[[sample_name]], sample_name)
  }

  message("Step 2-1 completed.") 
}

