require("RJSONIO")
require("paws")
require("zeallot")
require("ids")

reload_scdata_from_s3 <- function(pipeline_config, experiment_id) {
    s3 <- paws::s3(config=pipeline_config$aws_config)
    message(pipeline_config$source_bucket)
    message(paste(experiment_id, "r.rds", sep = "/"))
    
    c(body, ...rest) %<-% s3$get_object(
        Bucket = pipeline_config$source_bucket,
        Key = paste(experiment_id, "r.rds", sep = "/")
    )
    obj <- readRDS(rawConnection(body))
    return(obj)
}

send_output_to_api <- function(pipeline_config, input, plot_data_keys, output) {
    c(config, plot_data = plotData) %<-% output
    
    # upload output
    s3 <- paws::s3(config=pipeline_config$aws_config)
    id <- ids::uuid()
    output <- RJSONIO::toJSON(
        list(
            config = config,
            plotDataKeys = plot_data_keys,
            plotData = plot_data
        )
    )
    
    message("Uploading results to S3 bucket", pipeline_config$results_bucket, " at key ", id, "...")
    s3$put_object(
        Bucket = pipeline_config$results_bucket,
        Key = id,
        Body = charToRaw(output)
    )
    
    message("Sending to SNS topic ", pipeline_config$sns_topic)
    sns <- paws::sns(config=pipeline_config$aws_config)
    
    msg <- list(
        input = input,
        output = list(
            bucket = pipeline_config$results_bucket,
            key = id
        ),
        response = list(
            error = FALSE
        )
    )
    
    result <- sns$publish(
        Message = RJSONIO::toJSON(msg),
        TopicArn = pipeline_config$sns_topic,
        MessageAttributes = list(
            type = list(
                DataType = "String",
                StringValue = "PipelineResponse",
                BinaryValue = NULL
            )
        )
    )
    
    return(result$MessageId)
}

send_plot_data_to_s3 <- function(pipeline_config, experiment_id, output) { 
    plot_data <- output$plotData
    
    s3 <- paws::s3(config=pipeline_config$aws_config)
    
    plot_names <- names(plot_data)
    
    plot_data_keys = list()
    
    for(plot_data_name in names(plot_data)){
        id <- ids::uuid()
        
        plot_data_keys[[plot_data_name]] <- id
        
        output <- RJSONIO::toJSON(
            list(
                plotData = plot_data[[plot_data_name]]
            )
        )
        
        message("Uploading plotData to S3 bucket", pipeline_config$plot_data_bucket, " at key ", id, "...")
        s3$put_object(
            Bucket = pipeline_config$plot_data_bucket,
            Key = id,
            Body = charToRaw(output)
        )
    }
    
    return(plot_data_keys)
}

upload_matrix_to_s3 <- function(pipeline_config, experiment_id, data) {
    
    s3 <- paws::s3(config=pipeline_config$aws_config)
    
    object_key <- paste0(experiment_id, '/r.rds')
    
    count_matrix <- tempfile()
    saveRDS(data, file=count_matrix)
    
    message("Uploading updated count matrix to S3 bucket ", pipeline_config$processed_bucket, " at key ", object_key, "...")
    s3$put_object(
        Bucket = pipeline_config$processed_bucket,
        Key = object_key,
        Body = count_matrix
    )
    
    return(object_key)
}