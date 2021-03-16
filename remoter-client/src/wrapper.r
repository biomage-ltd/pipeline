require("RJSONIO")
require("paws")
require("zeallot")
require("ids")

load_config <- function() {
    config <- list(
        cluster_env = Sys.getenv("CLUSTER_ENV", "development"),
        sandbox_id = Sys.getenv("SANDBOX_ID", "default"),
        aws_account_id = Sys.getenv("AWS_ACCOUNT_ID", "242905224710"),
        aws_region = Sys.getenv("AWS_DEFAULT_REGION", "eu-west-1")
    )

    config[["aws_config"]] <- list(
        region = config$aws_region,
        credentials = list(
            creds = list(
                access_key_id = "mock-access-key",
                secret_access_key = "mock-secure-acces-key"
            )
        )
    )

    if(config$cluster_env == 'development') {
        config$aws_config[['endpoint']] <- 'http://host.docker.internal:4566'
        config$aws_account_id <- '000000000000'

        # This fixes a bug where paws would try to connect to InfraMock as if it was
        # a proper AWS endpoint, i.e. to http://bucket-name.host.docker.internal
        assignInNamespace("update_endpoint_for_s3_config",function(request) {
            return(request)
        },ns="paws.common")

    }

    config[["source_bucket"]] <- paste("biomage-source", config$cluster_env, sep = "-")
    config[["results_bucket"]] <- paste("worker-results", config$cluster_env, sep = "-")
    config[["sns_topic"]] <- paste("work-results", config$cluster_env, config$sandbox_id, sep = "-")
    config[["sns_topic"]] <- paste("arn:aws:sns", config$aws_region, config$aws_account_id, config$sns_topic, sep = ":")

    return(config)
}

pipeline_config <- load_config()

reload_from_s3 <- function(experiment_id) {
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


run_step <- function(task_name, scdata, config) {
    switch(task_name,
        test_fn = {
            import::from("/src/test_fn.r", task)
        },
        cellSizeDistribution = {
            import::from("/src/test_fn.r", task)
        },
        mitochondrialContent = {
            import::from("/src/test_fn.r", task)
        },
        classifier = {
            import::from("/src/test_fn.r", task)
        },
        numGenesVsNumUmis = {
            import::from("/src/test_fn.r", task)
        },
        doubletScores = {
            import::from("/src/doublet_scores.r", runDoubletScores)
        },
        dataIntegration = {
            import::from("/src/test_fn.r", task)
        },
        configureEmbedding = {
            import::from("/src/test_fn.r", task)
        },
        stop(paste("Invalid task name given:", task_name))
    )

    out <- task(scdata, config)
    return(out)
}

send_output_to_api <- function(input, output) {
    c(config, plot_data = plotData) %<-% output

    # upload output
    s3 <- paws::s3(config=pipeline_config$aws_config)
    id <- ids::uuid()
    output <- RJSONIO::toJSON(
        list(
            config = config,
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


wrapper <- function(input_json) {
    # Get data from state machine input.
    input <- RJSONIO::fromJSON(input_json)
    c(
        experiment_id = experimentId,
        task_name = taskName,
        config = config
    ) %<-% input

    if (!exists("scdata")) {
        message("No single-cell data has been loaded, reloading from S3...")

        # assign it to the global environment so we can
        # persist it across runs of the wrapper
        assign("scdata", reload_from_s3(experiment_id), pos = ".GlobalEnv")

        message("Single-cell data loaded.")
    }

    # call function to run and update global variable
    c(
        data, ...rest_of_results
    ) %<-% run_step(task_name, scdata, config)

    assign("scdata", data, pos = ".GlobalEnv")

    # send result to API
    message_id <- send_output_to_api(input, rest_of_results)

    return(message_id)
}

message("Wrapper loaded.")
