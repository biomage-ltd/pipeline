color_pool <- RJSONIO::fromJSON("/src/data-ingest/color_pool.json")

# set to local dirs for interactive dev
input_dir <- '/input'
output_dir <- '/output'

# creates the table information for samples
create_samples_table <- function(config, experiment_id, project_id) {
  # In samples_table we are going to add the core of the information
  samples_table <- list()

  # flag_filtered information
  df_prefiltered <- read.csv(file.path(output_dir, 'df_flag_filtered.txt'),
                             sep = '\t',
                             row.names = 'samples')

  samples <- row.names(df_prefiltered)
  samples_table$ids = paste0("sample-", samples)

  # For the current datasets it could happen that they are not in the gz format, so we leave the alternative tsv format.
  mime_options = c(
    "tsv" = "application/tsv",
    "gz" = "application/gzip",
    "mtx" = "application/mtx"
  )

  for (sample in samples) {

    prefiltered <- df_prefiltered[sample, 'flag_filtered'] == 'Filtered'

    # Identify datetime
    cdate <- mdate <- Sys.time()
    fnames <- list()

    # files that are not hidden
    sample_files <- file.path(
      sample,
      list.files(file.path(input_dir, sample))
    )

    # Iterate over each file to create the slot
    for (sample_file in sample_files) {

      fext <- tail(strsplit(sample_file, '[.]')[[1]], 1)

      fnames[[sample_file]] <- list(
        objectKey = '',
        name = sample_file,
        size = file.info(file.path(input_dir, sample_file))$size,
        mime = mime_options[[fext]],
        success = TRUE,
        error = FALSE
      )
    }

    # Add the whole information to each sample
    samples_table[[paste0("sample-", sample)]] <- list(
      name = sample,
      uuid = uuid::UUIDgenerate(),
      species = config$organism,
      type = config$input[['type']],
      createdDate = strftime(cdate, usetz = TRUE),
      lastModified = strftime(mdate, usetz = TRUE),
      complete = TRUE,
      error = FALSE,
      fileNames = sample_files,
      files = fnames,
      preFiltered = prefiltered
    )

  }


  return(list(
    "experimentId" = experiment_id,
    "samples" = samples_table,
    "projectUuid" = project_id
  ))
}


samples_sets <- function(){
  sample_annotations <- read.csv(file.path(output_dir, "samples-cells.csv"),
                                 sep = "\t",
                                 col.names = c("Cells_ID","Value"),
                                 na.strings = "None")

  cell_set <- list(key = "sample",
                   name = "Samples",
                   rootNode = TRUE,
                   children = list(),
                   type = "metadataCategorical")

  samples <- unique(sample_annotations$Value)

  for (sample in samples) {
    view <- sample_annotations[sample_annotations$Value == sample, "Cells_ID"]
    child <- list(key = paste0("sample-", sample),
                  name = sample,
                  color = color_pool[1],
                  cellIds = view)

    color_pool <- color_pool[-1]
    cell_set$children[[length(cell_set$children)+1]] <- child
  }

  return(cell_set)
}

# cell_sets fn for seurat metadata information
meta_sets <- function() {

  meta_annotations <- read.csv(file.path(output_dir, "metadata-cells.csv"), sep='\t')

  cell_set_list <- list()

  # The first column is the cells_id, the rest is the metadata information
  for (i in seq(2, ncol(meta_annotations))) {
    key <- name <- colnames(meta_annotations)[i]

    cell_set = list(
      "key" = key,
      "name" = name,
      "rootNode" = TRUE,
      "children" = list(),
      "type" = "metadataCategorical"
    )

    annot <- meta_annotations[[i]]

    for (value in unique(annot)) {
      view  <- meta_annotations[which(annot == value), 'cells_id']
      cell_set$children <- append(
        cell_set$children,
        list(
          "key" = paste(key, value, sep='-'),
          "name" = value,
          "color" = color_pool[1],
          "cellIds" = view)
      )

      color_pool <- color_pool[-1]
    }
    cell_set_list <- append(cell_set_list, cell_set)
  }
  return(cell_set_list)
}


task <- function(input, pipeline_config) {

  experiment_id <- input$experimentId
  project_id <- input$projectId
  sample_names <- input$sampleNames
  sample_uuids <- input$sampleUuids

  # save experiment_id for record-keeping
  writeLines(experiment_id, file.path(output_dir, "experiment_id.txt"))

  # read experiment config
  config <- RJSONIO::fromJSON(file.path(input_dir, "meta.json"))

  # read config related to QC pipeline
  config_dataProcessing <- RJSONIO::fromJSON(file.path(output_dir, "config_dataProcessing.json"))

  # Design cell_set scratchpad for DynamoDB
  scratchpad <- list(
    key = "scratchpad",
    name = "Scratchpad",
    rootNode = TRUE,
    children = list(),
    type = "cellSets"
  )

  # TODO: maybe we don't need samples_data
  samples_data <- create_samples_table(config, experiment_id, project_id)
  samples_set <- samples_sets()

  # Design cell_set meta_data for DynamoDB
  cell_sets <- list(scratchpad,samples_set)

  if ("metadata" %in% names(config))
    cell_sets <- append(cell_sets,meta_sets())
  
  cell_sets <- list(cellSets = cell_sets)

  print(paste("Experiment name is", config$name))

  experiment_data <- list(
    apiVersion = "2.0.0-data-ingest-seurat-rds-automated",
    experimentId = experiment_id,
    experimentName = config$name,
    meta = list(
      organism = config$organism,
      type = config$input[['type']]
    ),
    processingConfig = config_dataProcessing
  )

  # cell sets file to s3
  cell_sets_data <- RJSONIO::toJSON(cell_sets)

  put_object_in_s3(pipeline_config,
                   bucket = pipeline_config$cell_sets_bucket,
                   object = charToRaw(cell_sets_data),
                   key = experiment_id)

  # seurat object to s3
  put_object_in_s3(pipeline_config,
                   bucket = pipeline_config$source_bucket,
                   object = file.path(output_dir, 'experiment.rds'),
                   key = file.path(experiment_id, 'r.rds'))

  cluster_env <- pipeline_config$cluster_env
  print(sprintf("Experiment ID: %s uploaded to %s.", experiment_id, cluster_env))

  send_dynamodb_item_to_api(pipeline_config,
                            experiment_id = experiment_id,
                            table = pipeline_config$experiments_table,
                            item = experiment_data
                            task_name = "uploadToAWS")

  # samples data to dynamodb
  send_dynamodb_item_to_api(pipeline_config,
                            experiment_id = experiment_id,
                            table = pipeline_config$samples_table,
                            item = samples_data
                            task_name = "uploadToAWS")

  if (cluster_env == "production")
    print(sprintf("https://scp.biomage.net/experiments/%s/data-exploration", experiment_id))

}
