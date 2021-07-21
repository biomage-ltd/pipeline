download_cellranger <- function(input, pipeline_config) {
  project_id <- input$projectId
  sample_uuids <- input$sampleIds

  print_config(1, "Cellranger", input, pipeline_config, list())

  s3 <- paws::s3(config = pipeline_config$aws_config)

  fnames <- c("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")
  unlink("/input", recursive = TRUE)

  for (sample in sample_uuids) {
    for (fname in fnames) {
      gem_key <- file.path(project_id, sample, fname)

      message("GEM key")
      message(gem_key)

      # Preparing directories
      local_dir <- file.path("/input", sample)
      dir.create("/input")
      dir.create(local_dir)
      dir.create("/output")
      local_fpath <- file.path(local_dir, fname)

      # Download the file and store the output in a variable
      c(body, ...rest) %<-% s3$get_object(
        Bucket = pipeline_config$originals_bucket,
        Key = gem_key
      )

      # Write output to file
      writeBin(body, con = local_fpath)
    }
  }

  config <- list(
    name = input$experimentName,
    samples = input$sampleIds,
    organism = input$organism,
    input = list(type = input$input$type)
  )

  if ("metadata" %in% names(input)) {
    config$metadata <- input$metadata
  }

  exportJSON <- RJSONIO::toJSON(config)

  message("META file (config):")
  message(exportJSON)

  write(exportJSON, "/input/meta.json")

  message("Step 1 complete")
  return(list())
}
