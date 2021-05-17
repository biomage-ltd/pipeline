

# This function crate the table information for samples. As input it requires the experiment id and the config.
create_samples_table <- function(config, experimend_id) {
  # In samples_table we are going to add the core of the information
  samples_table <- list()

  # Getting flag_filtered information
  df_prefiltered <- read.csv('/output/df_flag_filtered.txt', sep = '\t', row.names = 'samples')
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
      list.files(file.path('/input', sample))
    )

    # Iterate over each file to create the slot
    for (sample_file in sample_files) {

      fext <- tail(strsplit(sample_file, '[.]')[[1]], 1)

      fnames[[sample_file]] <- list(
        objectKey = '',
        name = sample_file,
        size = file.info(file.path('/input',sample_file))$size,
        mime = mime_options[fext],
        success = TRUE,
        error = FALSE
      )
    }

    # Add the whole information to each sample
    samples_table[[paste0("sample-", sample)]] <- list(
      "name" = sample,
      "uuid" = uuid::UUIDgenerate(),
      "species" = config$organism,
      "type" = config$input$type,
      "createdDate" = strftime(cdate, usetz = TRUE),
      "lastModified" = strftime(mdate, usetz = TRUE),
      "complete" = TRUE,
      "error" = FALSE,
      "fileNames" = sample_files,
      "files" = fnames,
      "preFiltered" = prefiltered
    )



  }


  return(list(
    "experimentId" = experiment_id,
    "samples" = samples_table))

}


task <- function(experiment_id) {

  # save experiment_id for record-keeping
  writeLines(experiment_id, "/output/experiment_id.txt")


  config <- jsonlite::fromJSON("/input/meta.json")

  # read config related with QC pipeline
  config_dataProcessing <- jsonlite::fromJSON("/output/config_dataProcessing.json")

  # Design cell_set scratchpad for DynamoDB
  scratchpad = list(
    "key" = "scratchpad",
    "name" = "Scratchpad",
    "rootNode" = TRUE,
    "children" = c(),
    "type" = "cellSets"
  )

  samples_data = create_samples_table(config, experiment_id)
}
