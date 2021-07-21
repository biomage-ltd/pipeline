#  - Create a seurat object per sample
#  - Adding emtpyDrops,  scrublet and MT-content
#  - Create a file with flag_filtered

# LOADING SPARSE MATRIX AND CONFIGURATION
create_seurat <- function(input, pipeline_config) {
  message("reloading old matrices...")
  scdata_list <- readRDS("/output/pre-doublet-scdata_list.rds")

  message("Loading configuration...")
  config <- RJSONIO::fromJSON("/input/meta.json")

  print_config(5, "Create Seurat", input, pipeline_config, config)

  # Check which samples have been selected. Otherwiser we are going to use all of them.
  if (length(config$samples) > 0) {
    samples <- config$samples
  } else {
    samples <- names(scdata_list)
  }

  scdata_list <- scdata_list[samples]
  message("Creating Seurat Object...")

  sample_names <- unlist(input$sampleNames)
  names(sample_names) <- unlist(input$sampleIds)

  flag_filtered <- sapply(names(scdata_list), function(sample_id) {
    sample_name <- sample_names[sample_id]
    adding_metrics_and_annotation(scdata_list[[sample_id]], sample_id, sample_name, config)
  })

  df_flag_filtered <- data.frame(samples = samples, flag_filtered = ifelse(flag_filtered, "Filtered", "Unfiltered"))
  write.table(df_flag_filtered, "/output/df_flag_filtered.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

  message("Step 5 completed.")

  return(list())
}




# GETTING METADATA AND ANNOTATION


# adding_metrics_and_annotation function
#' @description We are going to process one sample in a seurat object. For it, we need to do:
#'  - Identify the possible metadata in the config file
#'  - Create a seurat object with the input of an sparse matrix
#'  - Annotate the gene in order to identify MT content in hsapiens and mmusculus
#'  - Computing MT-content
#'  - Getting scrublets
#'  - Getting emptyDrops
#'  - Identify flag_filtered
#'  - Save the rds with the seurat object for 1 sample
#' @param scdata Raw sparse matrix with the counts for one sample.
#' @param sample_name Name of the sample that we are preparing.
#' @param config Config of the project
#'
#' @return in the case that the input data was pre-filtered, we return a flag in order to disable the classifier filter.
#' This flag is going to be store in the dynamoDB inside the samples-table.
adding_metrics_and_annotation <- function(scdata, sample_id, sample_name, config, min.cells = 3, min.features = 10) {
  message("Converting into seurat object sample --> ", sample_id)

  metadata <- check_config(scdata, sample_id, sample_name, config)
  scdata <- Seurat::CreateSeuratObject(scdata, assay = "RNA", min.cells = min.cells, min.features = min.features, meta.data = metadata, project = config$name)

  organism <- config$organism

  # message("[", sample_id, "] \t finding genome annotations for genes...")
  # annotations <- gprofiler2::gconvert(
  # query = rownames(scdata), organism = organism, target="ENSG", mthreshold = Inf, filter_na = FALSE)

  annotations <- read.delim("/output/features_annotations.tsv")

  if (any(grepl("^mt-", annotations$name, ignore.case = T))) {
    message("[", sample_id, "] \t Adding MT information...")
    mt.features <- annotations$input[grep("^mt-", annotations$name, ignore.case = T)]
    mt.features <- mt.features[mt.features %in% rownames(scdata)]
    if (length(mt.features)) {
      scdata <- PercentageFeatureSet(scdata, features = mt.features, col.name = "percent.mt")
    }
  }

  if (is.null(scdata@meta.data$percent.mt)) scdata$percent.mt <- 0

  message("[", sample_id, "] \t Getting scrublet results...")
  scores <- get_doublet_score(sample_id)
  rownames(scores) <- scores$barcodes

  idt <- scores$barcodes[scores$barcodes %in% rownames(scdata@meta.data)]
  scdata@meta.data[idt, "doublet_scores"] <- scores[idt, "doublet_scores"]
  # Doublet class is the classifications that scDblFinder does to set the threshold of doublet_scores
  # (https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/2_scDblFinder.html#thresholding-and-local-calibration)
  scdata@meta.data[idt, "doublet_class"] <- scores[idt, "doublet_class"]

  message("[", sample_id, "] \t Adding emptyDrops...")
  file_ed <- paste("/output/pre-emptydrops-", sample_id, ".rds", sep = "")

  if (file.exists(file_ed)) {
    scdata@tools$flag_filtered <- FALSE
    message("\t \t getting emptyDrops results...")
    emptydrops_out <- readRDS(file = file_ed)

    emptydrops_out_df <- emptydrops_out %>%
      as.data.frame() %>%
      rlang::set_names(~ paste0("emptyDrops_", .)) %>%
      tibble::rownames_to_column("barcode")

    # adding emptydrops data to meta.data for later filtering (using left join)
    meta.data <- scdata@meta.data %>%
      tibble::rownames_to_column("barcode") %>%
      dplyr::left_join(emptydrops_out_df)
    rownames(meta.data) <- meta.data$barcode

    message("\t \t Adding emptyDrops scores information...")
    scdata@meta.data <- meta.data
    # previously (before joining into meta.data), results were just dumped as a additional slot
    # leaving the code here in case bugs arise from above solution
    # scdata@tools$CalculateEmptyDrops <- emptydrops_out
  } else {
    # Later on, when creating the config file, the enable will look the value of flag_filtered to   deactivate the classifier filter
    message("\t \t emptyDrops results not present, skipping...")
    scdata@meta.data$emptyDrops_FDR <- NA
    scdata@tools$flag_filtered <- TRUE
  }

  if (!dir.exists("/output/rds_samples")) {
    dir.create("/output/rds_samples")
  }

  message("[", sample_id, "] Saving R object...")
  saveRDS(scdata, file = paste("/output/rds_samples/", sample_id, ".rds", sep = ""), compress = FALSE)

  return(scdata@tools$flag_filtered)
}
