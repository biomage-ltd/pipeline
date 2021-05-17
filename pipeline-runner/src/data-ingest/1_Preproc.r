################################################
## 1_Preproc.r
##  - Read input folder {10x data}
##  - Prepare rds with a list of raw counts matrix per sample to compute emptyDrops 
##  - Prepare a csv per sample to run scrublet
################################################



suppressWarnings(library(Seurat))
suppressWarnings(library(magrittr))
suppressWarnings(library(parallel))
suppressWarnings(library(Matrix))
suppressWarnings(require(data.table))
suppressWarnings(library(RJSONIO))
suppressWarnings(library(MASS))


# checking_10x_structure
#' @description The design of the input data needs to be in a particular way. With this function we are going to check if
#' the input folder is designed correctly:
#' intput/
#' ------ sample_name/
#' ------------------- features.tsv.gz
#' ------------------- barcodes.tsv.gz
#' ------------------- matrix.mtx.gz
#'
#'cell ranger output
#'https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
#'STAR solo conventions (drop in replacement for cell ranger):
#'https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
#' @param samples samples to check
#'
#' @return TRUE if the design is correct FALSE otherwise
check_10x_input <- function(samples){
  
  # Cell ranger v2 does in fact produce non gunzipped outputs, but since during
  # data-upload gzip compresses all non-compressed file, so we expect files with *.gz ending here
  cell_ranger_v2 <- c("genes.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")
  cell_ranger_v3 <- c("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")

  for(sample in samples){
    check_v2 <- all(cell_ranger_v2%in%list.files(paste("/input", sample, sep="/"), full.names = F))
    check_v3 <- all(cell_ranger_v3%in%list.files(paste("/input", sample, sep="/"), full.names = F))

    if(!check_v2 & !check_v3){
      return(FALSE)
    }
  }

  if(check_v2)
    message("Version of Cell Ranger: V2")

  if(check_v3)
    message("Version of Cell Ranger: V3")

  return(TRUE)

}

# Read10X only support CellRanger V2 for ungz version of the file. However, we will
# store all the files with the gz version. This function add/remove the gz from given files.

# create_dataframe function
#' @description read matrix based on config. Possibles input
#'      - 10x scdata
#'      - table
#' @param config experiment settings.
#'
#' @return list with an element per sample with dgCMatrix with genes as rows and cells as column

create_dataframe <- function(config){
  data_type <- config$input["type"]
  scdata <- list()

  # Check config format. It requires to have a 'samples' key. For the multisample experiment, it needs to be the same
  # as the name of the folders that are inside the input folder. In the case of unisample it should be [].
  if(!"samples"%in%names(config))
    stop("The format of the config is wrong. There should be a category for the sample.
     Current content of config:", names(config))
  
  # If no samples have been added, we will take all the samples of the project. 
  if(length(config$samples)==0){
    warning("All samples in input folder will be included in the analysis!")
    samples <- list.dirs("/input", full.names = F)[-1]
  }else{
    samples <- config$samples

    if (!all(samples%in%list.dirs("/input", full.names = F)))
      stop("Check samples to be used in the analysis, 
      since there are some of them that hasn't got a folder with the files: ",
           samples[!samples%in%list.dirs("/input", full.names = F)])
  }
  
  message("Samples to include in the analysis: ", paste(samples, collapse = " - "))
  
  if (data_type == "10x"){
    message("Loading 10x data set from input folder.")

    # Checking design
    if(!check_10x_input(samples)){
      stop("Please! Check in files inside the input folder.", 
      "There should be the files {genes.tsv, matrix.mtx, barcodes.tsv} or {features.tsv.gz, barcodes.tsv.gz and matrix.mtx.gz}")
    }

    annotation_features <- list()
    # overall feature annotation is derived from input data saved in genes.tsv features.tsv.gz
    # since each sample only carries a subset of annotation for it's expressed genes, the annotation for all samples is merged.
    # this is an excerpt of features.tsv.gz
    # ENSG00000237613 | FAM138A | Gene Expression
    # ENSG00000186092 | OR4F5 | Gene Expression
    # More information about genes.tsv features.tsv.gz: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices

    for(sample in samples){
      sample_dir <- file.path('/input', sample)


      sample_fpaths <-  list.files(sample_dir, full.names = TRUE)
      message("Reading files --> ", sample_fpaths)

      # If we are in CellRanger V2 we need to remove gz to the files in order to use the function Read10X
      is_v2 <- any(grepl("genes.tsv(.gz)?$", sample_fpaths))
      if (is_v2) {
          new_fpaths <- gsub(".gz$", "", sample_fpaths)
          file.rename(sample_fpaths, new_fpaths)
          sample_fpaths <- new_fpaths
      }

      scdata[[sample]] <- Read10X(sample_dir, gene.column = 1)  
      annot_fpath <- sample_fpaths[grepl('genes.tsv$|features.tsv.gz$', sample_fpaths)] 
      annotation_features[[sample]] <- read.delim(annot_fpath, header = FALSE)

      message(
        paste(
          "Found", nrow(scdata[[sample]]), "genes and", ncol(scdata[[sample]]), "cells in sample", sample, "."
        )
      )

      # Now, we can add again the gz suffix
      if (is_v2) file.rename(sample_fpaths, paste0(sample_fpaths, '.gz'))

    }
    annotation_features_df <- unique(do.call('rbind', annotation_features))
    annotation_features_df <- annotation_features_df[, c(1, 2)]
    colnames(annotation_features_df) <- c("input", "name")
    write.table(annotation_features_df, "/output/features_annotations.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
  
  # So far, we haven't tested input format in table type. 
  if (data_type == "table") {
    message(paste("Loading table-type data set from", path))
    for(sample in samples){
      scdata[[sample]] <- as.matrix(read.table(paste("/input", sample, sep = "/")))
      message(
        paste(
          "Found", nrow(scdata$raw[[sample]]), "genes and", ncol(scdata$raw[[sample]]), "cells in sample", sample, "."
        )
      )
    }
  }
  
  return(scdata)
}

task <- function(input,pipeline_config) {
  message("Loading configuration...")
  config <- RJSONIO::fromJSON("/input/meta.json")

  print("Config:")
  print(config)
  print("Current working directory:")
  print(getwd())
  print("Experiment folder status:")
  print(list.files(paste("/input",sep = "/"),all.files=TRUE,full.names=TRUE,recursive=TRUE))


  # We include in variable scdata_list all the sparse matrix per sample
  message("Creating raw dataframe...")
  scdata_list <- create_dataframe(config)

  # We store the raw scdata_list since for the emptyDrops since to compute the background we cannot remove any cells. 
  message("Exporting raw scdata for emptyDrops...")
  saveRDS(scdata_list, file = "/output/pre-doublet-scdata_list.rds", compress = FALSE)

  message("Step 1 completed.")
  print(list.files(paste("/output",sep = "/"),all.files=TRUE,full.names=TRUE,recursive=TRUE))
}
