require("RJSONIO")
require("paws")
require("zeallot")
require("ids")

task <- function(input,pipeline_config) {
    # you are receiving the full sample object instead of sample names
    project_id <- input$projectId
    sample_names <- input$sampleNames # extract sample names from samples object
    sample_uuids <- input$sampleUuids
    message("download2")
    s3 <- paws::s3(config=pipeline_config$aws_config)
    message(pipeline_config$originals_bucket)
    
    fnames <- c('features.tsv.gz', 'barcodes.tsv.gz', 'matrix.mtx.gz')
    unlink("/input",recursive=TRUE)
    for (sample in sample_uuids) {
        for (fname in fnames) {
            gem_key <- file.path(project_id, sample, fname)
            message(gem_key)
            sample_name = sample_names[[match(sample,sample_uuids)]]
            #Preparing directories
            local_dir <- file.path('/input',sample_name)
            #unlink(local_dir, recursive = TRUE)
            dir.create('/input')
            dir.create(local_dir)
            dir.create("/output")
            local_fpath <- file.path(local_dir,fname)
            
            message("bucket")
            message( pipeline_config$originals_bucket)
            message("file")
            message(gem_key)
            # Download the file and store the output in a variable
            c(body, ...rest) %<-% s3$get_object(
                #Bucket = pipeline_config$originals_bucket,
                Bucket = pipeline_config$originals_bucket,
                Key = gem_key
            )
            
            # Write output to file
            writeBin(body, con = local_fpath)
        }
    }
    #ownload meta.json
    meta_key = file.path(project_id, "meta.json")
    message(paste("File: ",meta_key))
    c(body, ...rest) %<-% s3$get_object(
    Bucket = pipeline_config$originals_bucket,
    
    Key = meta_key
    )
    writeBin(body, con = "/input/meta.json")
}