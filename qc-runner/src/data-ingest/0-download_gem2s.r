require("RJSONIO")
require("paws")
require("zeallot")
require("ids")

task <- function(input,pipeline_config) {
    project_id <- input$projectId
    sample_names <- input$sampleNames

    s3 <- paws::s3(config=pipeline_config$aws_config)
    message(pipeline_config$originals_bucket)
    
    fnames <- c('features.tsv.gz', 'barcodes.tsv.gz', 'matrix.mtx.gz')
    
    for (sample in sample_names) {
        for (fname in fnames) {
            gem_key <- file.path(project_id, sample, fname)
            message(gem_key)
            
            #Preparing directories
            local_dir <- file.path('/input',project_id,sample)
            unlink(local_dir, recursive = TRUE)
            dir.create(local_dir)
            dir.create("/output")
            local_fpath <- file.path(local_dir,fname)
            
            
            # Download the file and store the output in a variable
            c(body, ...rest) %<-% s3$get_object(
                Bucket = pipeline_config$originals_bucket,
                Key = gem_key
            )
            
            # Write output to file
            writeBin(body, con = local_fpath)
        }
    }
    #Download meta.json
    c(body, ...rest) %<-% s3$get_object(
    Bucket = pipeline_config$originals_bucket,
    Key = "meta.json"
    )
    writeBin(body, con = "/input/meta.json")
}