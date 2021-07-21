#  - Merging the samples for the current experiment
#  - Adding metadata: cellsId, color_pool and gene annotation
#  - Preparing dataProcessing json file


prepare_experiment <- function(input, pipeline_config) {
  message("Loading configuration...")
  config <- RJSONIO::fromJSON("/input/meta.json")

  print_config(6, "Prepare experiment", input, pipeline_config, config)

  # Check which samples have been selected. Otherwiser we are going to use all of them.
  if (length(config$samples) > 0) {
    samples <- config$samples
  } else {
    samples <- gsub("\\..*", "", list.files("/output/rds_samples"))
  }

  message("Reloading samples rds for current experiment...")
  scdata_list <- list()
  for (sample in samples) {
    scdata_list[[sample]] <- readRDS(paste("/output/rds_samples/", sample, ".rds", sep = ""))
  }

  # Merging samples and adding a prefix with the sample name. In pipeline we grep in barcodes to filter by sample.
  if (length(scdata_list) == 1) {
    scdata <- scdata_list[[1]]
    scdata <- RenameCells(object = scdata, add.cell.id = names(scdata_list)[1])
  } else {
    scdata <- merge(scdata_list[[1]], y = scdata_list[-1], add.cell.ids = c(samples))
  }

  message("Storing gene annotations...")
  organism <- config$organism
  annotations <- read.delim("/output/features_annotations.tsv")

  # In order to avoid duplicated genes names, we are going to add the ENSEMBL ID for those
  # genes that are duplicated (geneNameDuplicated-ENSEMBL)
  gname <- annotations$name
  # Keep original name in 'original_name' variable
  annotations$original_name <- gname
  is.dup <- duplicated(gname) | duplicated(gname, fromLast = TRUE)
  annotations$name[is.dup] <- paste(gname[is.dup], annotations$input[is.dup], sep = " - ")

  # Ensure index by rownames in scdata
  annotations <- annotations[match(rownames(scdata), annotations$input), ]
  rownames(annotations) <- annotations$input

  scdata@misc[["gene_annotations"]] <- annotations

  message("Storing cells id...")
  # Keeping old version of ids starting from 0
  scdata$cells_id <- 0:(nrow(scdata@meta.data) - 1)

  message("Storing color pool...")
  # We store the color pool in a slot in order to be able to access it during configureEmbedding
  scdata@misc[["color_pool"]] <- get_color_pool()
  message("Stored pool")

  scdata@misc[["experimentId"]] <- input$experimentId
  scdata@misc[["ingestionDate"]] <- Sys.time()


  if ("metadata" %in% names(config)) {
    message("Storing metadata...")
    scdata@misc[['meta_vars']] <- names(config$metadata)
  }

  # CHECK FILTERED DATA
  # -

  df_flag_filtered <- read.delim("/output/df_flag_filtered.txt")
  any_filtered <- "Filtered" %in% df_flag_filtered$flag_filtered
  message("saved filtered flag")

  # TEST OBJECT
  # -

  # test_object(scdata)

  # SAVING FILES
  # -

  message("saving R object...")
  saveRDS(scdata, file = "/output/experiment.rds", compress = FALSE)

  write.table(
    colnames(scdata),
    file = "/output/r-out-cells.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
  )

  write.table(
    scdata@misc[["gene_annotations"]][scdata@misc[["gene_annotations"]]$input %in% rownames(scdata), ],
    file = "/output/r-out-annotations.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
  )

  print(scdata)

  # DATA PROCESSING
  # --
  # [HARDCODED]
  config.cellSizeDistribution <- list(
    enabled = "true",
    auto = "true",
    filterSettings = list(minCellSize = 1080, binStep = 200)
  )

  config.mitochondrialContent <- list(
    enabled = "true", auto = "true",
    filterSettings = list(
      method = "absolute_threshold",
      methodSettings = list(absolute_threshold = list(maxFraction = 0.1, binStep = 0.05))
    )
  )

  config.classifier <- list(
    enabled = tolower(as.character(!any_filtered)), # emptyDrops results not present
    auto = "true",
    filterSettings = list(FDR = 0.01)
  )

  config.numGenesVsNumUmis <- list(
    enabled = "true",
    auto = "true",
    filterSettings = list(
      regressionType = "gam",
      regressionTypeSettings = list("gam" = list(p.level = 0.001))
    )
  )

  config.doubletScores <- list(
    enabled = "true",
    auto = "true",
    filterSettings = list(probabilityThreshold = 0.5, binStep = 0.05)
  )

  # BE CAREFUL! The method is based on config.json. For multisample only seuratv4, for unisample LogNormalize
  # hardcoded because unisample check is performed in dataIntegration
  identified.method <- "harmony"
  config.dataIntegration <- list(
    dataIntegration = list(
      method = identified.method,
      methodSettings = list(
        seuratv4 = list(numGenes = 2000, normalisation = "logNormalize"),
        unisample = list(numGenes = 2000, normalisation = "logNormalize"),
        harmony = list(numGenes = 2000, normalisation = "logNormalize"),
        fastmnn = list(numGenes = 2000, normalisation = "logNormalize")
      )
    ),
    dimensionalityReduction = list(method = "rpca", numPCs = 30, excludeGeneCategories = c())
  )

  config.configureEmbedding <- list(
    embeddingSettings = list(
      method = "umap",
      methodSettings = list(
        umap = list(minimumDistance = 0.3, distanceMetric = "euclidean"),
        tsne = list(
          perplexity = min(30, ncol(scdata) / 100),
          learningRate = max(200, ncol(scdata) / 12)
        )
      )
    ),
    clusteringSettings = list(
      method = "louvain",
      methodSettings = list(louvain = list(resolution = 0.8))
    )
  )

  mitochondrial_config_to_duplicate <- list(
    auto = "true",
    filterSettings = list(
      method = "absolute_threshold",
      methodSettings = list(absolute_threshold = list(maxFraction = 0.1, binStep = 0.05))
    )
  )

  classifier_config_to_duplicate <- list(
    enabled = tolower(as.character(!any_filtered)), # emptyDrops results not present
    auto = "true",
    filterSettings = list(FDR = 0.01)
  )

  samples <- scdata$samples

  # Add a copy of the same base config for each sample for classifier and mitochondrialContent
  config.classifier <- duplicate_config_per_sample(classifier_config_to_duplicate, config.classifier, samples)
  config.mitochondrialContent <- duplicate_config_per_sample(mitochondrial_config_to_duplicate, config.mitochondrialContent, samples)

  # Compute for multisample and unisample
  config.cellSizeDistribution <- add_custom_config_per_sample(cellSizeDistribution_config, config.cellSizeDistribution, scdata)
  config.numGenesVsNumUmis <- add_custom_config_per_sample(numGenesVsNumUmis_config, config.numGenesVsNumUmis, scdata)
  config.doubletScores <- add_custom_config_per_sample(doubletScores_config, config.doubletScores, scdata)

  # When we remove the steps from data-ingest we need to change here the default config.
  # Save config for all steps.
  config <- list(
    cellSizeDistribution = config.cellSizeDistribution,
    mitochondrialContent = config.mitochondrialContent,
    classifier = config.classifier,
    numGenesVsNumUmis = config.numGenesVsNumUmis,
    doubletScores = config.doubletScores,
    dataIntegration = config.dataIntegration,
    configureEmbedding = config.configureEmbedding
  )


  # Export to json
  exportJson <- RJSONIO::toJSON(config, pretty = TRUE)
  # The RJSONIO library add '' to boolean keys, so we will remove them.
  exportJson <- gsub('\"true\"', "true", exportJson)
  exportJson <- gsub('\"false\"', "false", exportJson)
  # Tranform null into []
  exportJson <- gsub("null", "[]", exportJson)
  message("config file...")
  write(exportJson, "/output/config_dataProcessing.json")

  message("Step 6 completed.")

  return(list())
}


cellSizeDistribution_config <- function(scdata, config) {
  minCellSize <- generate_default_values_cellSizeDistribution(scdata, config, 1e2)
  config$filterSettings$minCellSize <- minCellSize
  return(config)
}

# To identify intelligently the treshold we are going to use the logic inside scDblFinder, which creates a classification
# (singlet our doublet) [ref: https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/2_scDblFinder.html#thresholding-and-local-calibration]
# To set the auto value we are going to use as a threshold the maximun score that is given to a singlet.

doubletScores_config <- function(scdata, config) {
  # Maximum score that has a singlet
  probabilityThreshold <- max(scdata$doublet_scores[scdata$doublet_class == "singlet"], na.rm = TRUE)
  config$filterSettings$probabilityThreshold <- probabilityThreshold

  return(config)
}


numGenesVsNumUmis_config <- function(scdata, config) {
  # Sensible values are based on the funciton "gene.vs.molecule.cell.filter" from the pagoda2 package
  p.level <- min(0.001, 1 / ncol(scdata))
  config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]]$p.level <- p.level

  return(config)
}


# SAVING CONFIG FILE
#
# We are going to store the final config to config_dataProcessing.json, in order to upload to dynamoDB.
# The unisample experiments does not require any change, but for the multisample experiment we need
# to add the filtering parameter for each sample (only in the steps that is required.)
# We are going to differentiate in samples only in the steps:
# --> cellSizeDistribution
# --> numGenesVsNumUmis
# --> doubletScores
#
# For both of them, we will run again the step fn for each sample (samples names are stored in metadata type)

# Function to recompute the step fn and store the new config of each sample inside the latest config file
# We need to iterate per sample and compute separately the step fn.
# Example of structure:
# {
# “filterSettings”: {
#   “probabilityThreshold”: 0.2,
#   “binStep”: 0.05
# },
# “sample-KO”: {
#   “filterSettings”: {
#     “probabilityThreshold”: 0.1,
#     “binStep”: 100
#   }
# },
# “sample-WT1": {
#                     “filterSettings”: {
#                         “probabilityThreshold”: 0.1,
#                         “binStep”: 45
#                     }
#                 }
# }

duplicate_config_per_sample <- function(step_config, config, samples) {
  for (sample in unique(samples)) {
    config[[sample]] <- step_config
    config[[sample]]$defaultFilterSettings <- step_config$filterSettings
  }

  return(config)
}

add_custom_config_per_sample <- function(step_fn, config, scdata) {

  # We upadte the config file, so to be able to access the raw config we create a copy
  config.raw <- config

  samples <- scdata$samples

  for (sample in unique(samples)) {
    # Downsample the seurat object to a unisample experiment
    scdata_sample <- scdata[, samples %in% sample]
    # Run the step fun with the unisample experiment and keep the config result
    result_config <- step_fn(scdata_sample, config.raw)
    # Update config with the unisample thresholds
    config[[sample]] <- result_config

    # Add auto settings
    config[[sample]]$defaultFilterSettings <- result_config$filterSettings
  }

  return(config)
}
