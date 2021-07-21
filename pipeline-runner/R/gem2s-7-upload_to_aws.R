upload_to_aws <- function(input, pipeline_config) {
  # set to local dirs for interactive dev
  input_dir <- "/input"
  output_dir <- "/output"

  experiment_id <- input$experimentId
  project_id <- input$projectId

  # save experiment_id for record-keeping
  writeLines(experiment_id, file.path(output_dir, "experiment_id.txt"))

  # read experiment config
  config <- RJSONIO::fromJSON(file.path(input_dir, "/meta.json"))

  print_config(7, "Upload to AWS", input, pipeline_config, config)

  # read config related to QC pipeline
  config_dataProcessing <- RJSONIO::fromJSON(file.path(output_dir, "config_dataProcessing.json"))


  experiment_data <- list(
    apiVersion = "2.0.0-data-ingest-seurat-rds-automated",
    experimentId = experiment_id,
    experimentName = config$name,
    meta = list(
      organism = config$organism,
      type = config$input[["type"]]
    ),
    processingConfig = config_dataProcessing
  )

  # need to create initial cell sets to preview plots
  scdata <- readRDS(file.path(output_dir, "experiment.rds"))
  update_cell_sets(scdata, experiment_id, pipeline_config)

  # seurat object to s3
  put_object_in_s3_multipart(pipeline_config,
    bucket = pipeline_config$source_bucket,
    object = file.path(output_dir, "experiment.rds"),
    key = file.path(experiment_id, "r.rds")
  )

  cluster_env <- pipeline_config$cluster_env
  print(sprintf("Experiment ID: %s uploaded to %s.", experiment_id, cluster_env))

  data <- list(
    item = experiment_data,
    table = pipeline_config$experiments_table
  )

  if (cluster_env == "production") {
    print(sprintf("https://scp.biomage.net/experiments/%s/data-exploration", experiment_id))
  }

  return(data)
}

