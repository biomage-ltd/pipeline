require("RJSONIO")
require("paws")
require("zeallot")
require("ids")

# if not attached can cause errors when accessing metadata
require("Seurat")

load_config <- function(development_aws_server) {
    config <- list(
        cluster_env = Sys.getenv("CLUSTER_ENV", "development"),
        sandbox_id = Sys.getenv("SANDBOX_ID", "default"),
        aws_account_id = Sys.getenv("AWS_ACCOUNT_ID", "242905224710"),
        aws_region = Sys.getenv("AWS_DEFAULT_REGION", "eu-west-1"),
        pod_name = Sys.getenv("K8S_POD_NAME", "local"),
        activity_arn = Sys.getenv("ACTIVITY_ARN", "")
    )

    config[["aws_config"]] <- list(
        region = config$aws_region
    )

    if(config$cluster_env == 'development') {
        config$aws_config[['endpoint']] <- sprintf("http://%s:4566", development_aws_server) # DOCKER_GATEWAY_HOST
        config$aws_config[['credentials']] <- list(
            creds = list(
                access_key_id = "mock-access-key",
                secret_access_key = "mock-secure-acces-key"
            )
        )
        config$aws_account_id <- '000000000000'

        # This fixes a bug where paws would try to connect to InfraMock as if it was
        # a proper AWS endpoint, i.e. to http://bucket-name.host.docker.internal
        assignInNamespace("update_endpoint_for_s3_config",function(request) {
            return(request)
        },ns="paws.common")

    }

    config[["source_bucket"]] <- paste("biomage-source", config$cluster_env, sep = "-")
    config[["processed_bucket"]] <- paste("processed-matrix", config$cluster_env, sep = "-")
    config[["results_bucket"]] <- paste("worker-results", config$cluster_env, sep = "-")
    config[["plot_data_bucket"]] <- paste("plots-tables", config$cluster_env, sep = "-")
    config[["sns_topic"]] <- paste("work-results", config$cluster_env, config$sandbox_id, sep = "-")
    config[["sns_topic"]] <- paste("arn:aws:sns", config$aws_region, config$aws_account_id, config$sns_topic, sep = ":")

    return(config)
}

reload_from_s3 <- function(pipeline_config, experiment_id) {
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

run_step <- function(scdata, config, task_name, sample_id) {
    switch(task_name,
        cellSizeDistribution = {
            import::here("/src/cellSizeDistribution.r", task)
        },
        mitochondrialContent = {
            import::here("/src/mitochondrialContent.r", task)
        },
        classifier = {
            import::here("/src/classifier.r", task)
        },
        numGenesVsNumUmis = {
            import::here("/src/numGenesVsNumUmis.r", task)
        },
        doubletScores = {
            import::here("/src/doubletScores.r", task)
        },
        dataIntegration = {
            import::here("/src/dataIntegration.r", task)
        },
        configureEmbedding = {
            import::here("/src/configureEmbedding.r", task)
        },
        stop(paste("Invalid task name given:", task_name))
    )
    out <- task(scdata, config, task_name, sample_id)
    return(out)
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


wrapper <- function(input_json) {
    # Get data from state machine input.
    input <- RJSONIO::fromJSON(input_json, simplify = FALSE)
    
    str(input)

    experiment_id = input$experimentId
    task_name = input$taskName
    config = input$config
    server = input$server
    upload_count_matrix = input$uploadCountMatrix
    sample_id = input$sampleUuid

    if (sample_id != "") {
        config <- config[[sample_id]]
        input$config <- config
    }
    
    input <- input[names(input) != "server"]

    pipeline_config <- load_config(server)

    if (!exists("scdata")) {
        message("No single-cell data has been loaded, reloading from S3...")

        # assign it to the global environment so we can
        # persist it across runs of the wrapper
        assign("scdata", reload_from_s3(pipeline_config, experiment_id), pos = ".GlobalEnv")

        message("Single-cell data loaded.")
    }

    # call function to run and update global variable
    c(
        data, ...rest_of_results
    ) %<-% run_step(scdata, config, task_name, sample_id)

    assign("scdata", data, pos = ".GlobalEnv")

    # upload plot data result to S3
    plot_data_keys <- send_plot_data_to_s3(pipeline_config, experiment_id, rest_of_results)

    # Uplaod count matrix data
    if(upload_count_matrix) {
        object_key <- upload_matrix_to_s3(pipeline_config, experiment_id, scdata)
        message('Count matrix uploaded to ', pipeline_config$processed_bucket, ' with key ',object_key)
    }

    # send result to API
    message_id <- send_output_to_api(pipeline_config, input, plot_data_keys, rest_of_results)

    return(message_id)
}

init <- function() {
    pipeline_config <- load_config('host.docker.internal')
    states <- paws::sfn(config=pipeline_config$aws_config)

    repeat {
        c(taskToken, input) %<-% states$get_activity_task(
            activityArn = pipeline_config$activity_arn,
            workerName = pipeline_config$pod_name
        )

        if(is.null(taskToken) || !length(taskToken) || taskToken == "") {
            message('No input received during last poll, shutting down...')
            quit('no')
        }

        tryCatch(
            expr = {
                message('Input ', input, ' found')        
                wrapper(input)

                states$send_task_success(
                    taskToken = taskToken,
                    output = '{}'
                )
            },
            error = function(e){ 
                message("Error: ",e)
                states$send_task_failure(
                    taskToken = taskToken,
                    error = 'R error',
                    cause = e
                )
            }
        )
    }
}

init()
