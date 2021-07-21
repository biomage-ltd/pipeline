# if Seurat not attached can cause errors when accessing metadata
library(Seurat)
library(zeallot)

# increase maxSize from the default of 500MB to 32GB
options(future.globals.maxSize = 32 * 1024 * 1024^2)

for (f in list.files('R', '.R$', full.names = TRUE)) source(f)
library(magrittr)

load_config <- function(development_aws_server) {
    config <- list(
        cluster_env = Sys.getenv("CLUSTER_ENV", "development"),
        sandbox_id = Sys.getenv("SANDBOX_ID", "default"),
        aws_account_id = Sys.getenv("AWS_ACCOUNT_ID", "242905224710"),
        aws_region = Sys.getenv("AWS_DEFAULT_REGION", "eu-west-1"),
        pod_name = Sys.getenv("K8S_POD_NAME", "local"),
        activity_arn = Sys.getenv("ACTIVITY_ARN", ""),
        debug_config = list(
            step = Sys.getenv("DEBUG_STEP", ""),
            path = Sys.getenv("DEBUG_PATH", "")
        )
    )

    config[["aws_config"]] <- list(region = config$aws_region)

    # running in linux needs the IP of the host to work. If it is set as an environment variable (by makefile) honor it instead of the
    # provided parameter
    overriden_server = Sys.getenv("HOST_IP", "")
    if (overriden_server != "") {
        development_aws_server = overriden_server
    }

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
        assignInNamespace("update_endpoint_for_s3_config", function(request) {
            return(request)
        }, ns="paws.common")

    }

    config[["samples_table"]] <- paste("samples", config$cluster_env, sep = "-")
    config[["experiments_table"]] <- paste("experiments", config$cluster_env, sep = "-")
    config[["originals_bucket"]] <- paste("biomage-originals", config$cluster_env, sep = "-")
    config[["source_bucket"]] <- paste("biomage-source", config$cluster_env, sep = "-")
    config[["processed_bucket"]] <- paste("processed-matrix", config$cluster_env, sep = "-")
    config[["results_bucket"]] <- paste("worker-results", config$cluster_env, sep = "-")
    config[["plot_data_bucket"]] <- paste("plots-tables", config$cluster_env, sep = "-")
    config[["cell_sets_bucket"]] <- paste("cell-sets", config$cluster_env, sep = "-")
    config[["sns_topic"]] <- paste("work-results", config$cluster_env, config$sandbox_id, sep = "-")
    config[["sns_topic"]] <- paste("arn:aws:sns", config$aws_region, config$aws_account_id, config$sns_topic, sep = ":")


    return(config)
}

run_processing_step <- function(scdata, config, task_name, sample_id, debug_config) {

    # vector of task functions named by task name
    tasks <- list(
        'classifier' = filter_emptydrops,
        'cellSizeDistribution' = filter_low_cellsize,
        'mitochondrialContent' = filter_high_mito,
        'numGenesVsNumUmis' = filter_gene_umi_outlier,
        'doubletScores' = filter_doublets,
        'dataIntegration' = integrate_scdata,
        'configureEmbedding' = embed_and_cluster
    )

    if (!task_name %in% names(tasks)) stop('Invalid task: ', task_name)

    handle_debug(scdata, config, task_name, sample_id, debug_config)

    # print info
    task <- tasks[[task_name]]
    message("Running: ", task_name)
    message("Config:")
    str(config)

    # run task and time it
    tstart <- Sys.time()
    out <- task(scdata, config, sample_id, task_name)
    ttask <- format(Sys.time()-tstart, digits = 2)
    message("⏱️ Time to complete ", task_name, " for sample ", sample_id, ": ", ttask)

    # filter specific info
    is.filter <- task_name %in% names(tasks)[1:5]
    if (is.filter) {
        message("Cells per sample before filter for sample: ", sample_id)
        print(table(scdata$samples))

        message("Cells per sample after filter for sample: ", sample_id)
        print(table(out$data$samples))
    }

    return(out)
}

#
# input needs:
# sampleIds for all samples
# sampleNames for all samples
# projectId
# pipeline_config as defined by load_config
#

run_gem2s_step <- function(task_name, input, pipeline_config) {

    # list of task functions named by task name
    tasks <- list(
        'downloadGem' = download_cellranger,
        'preproc' = load_cellranger,
        'emptyDrops' = run_emptydrops,
        'doubletScores' = score_doublets,
        'createSeurat' = create_seurat,
        'prepareExperiment' = prepare_experiment,
        'uploadToAWS' = upload_to_aws
    )

    if (!task_name %in% names(tasks)) stop("Invalid task name given: ", task_name)
    message("Starting task: ", task_name)
    task <- tasks[[task_name]]

    tstart <- Sys.time()
    data <- task(input, pipeline_config)
    ttask <- format(Sys.time()-tstart, digits = 2)
    message("⏱️ Time to complete ", task_name, ": ", ttask)

    return(data)
}

call_gem2s <- function(task_name, input, pipeline_config) {
    experiment_id <- input$experimentId

    data %<-% run_gem2s_step(task_name, input, pipeline_config)
    message_id <- send_gem2s_update_to_api(pipeline_config, experiment_id, task_name, data)

    return(message_id)
}

#
# call_data_processing
# Runs step @task_name of the data processing pipeline, send plot data to s3 and output message to api.
# IN task_name: str, name of the step.
# IN input: json str, config settings for all samples, current sample uuid, uploadCountMatrix (whether or not to upload after step)
# IN pipeline_config: json str, environment config resulting from the load_config function
#
call_data_processing <- function(task_name, input, pipeline_config) {

    experiment_id <- input$experimentId
    config <- input$config
    upload_count_matrix <- input$uploadCountMatrix
    upload_cell_sets <- input$uploadCellSets
    sample_id <- input$sampleUuid
    debug_config <- pipeline_config$debug_config

    if (sample_id != "") {
        config <- config[[sample_id]]
        input$config <- config
    }


    if (!exists("scdata")) {
        message("No single-cell data has been loaded, reloading from S3...")

        # assign it to the global environment so we can
        # persist it across runs of the wrapper
        assign("scdata", reload_scdata_from_s3(pipeline_config, experiment_id), pos = ".GlobalEnv")

        message("Single-cell data loaded.")
    }

    # call function to run and update global variable
    c(
        data, ...rest_of_results
    ) %<-% run_processing_step(scdata, config, task_name, sample_id, debug_config)

    assign("scdata", data, pos = ".GlobalEnv")

    # upload plot data result to S3
    plot_data_keys <- send_plot_data_to_s3(pipeline_config, experiment_id, rest_of_results)

    # update and upload cellSets
    if (upload_cell_sets) {
        message('Uploading cell sets to S3')
        update_cell_sets(scdata, experiment_id, pipeline_config)
    }

    # Uplaod count matrix data
    if(upload_count_matrix) {
        object_key <- upload_matrix_to_s3(pipeline_config, experiment_id, scdata)
        message('Count matrix uploaded to ', pipeline_config$processed_bucket, ' with key ',object_key)
    }

    # send result to API
    message_id <- send_output_to_api(pipeline_config, input, plot_data_keys, rest_of_results)

    return(message_id)
}

#
# Wrapper(input_json)
# IN input_json: json input from message. Input should have:
# taskname, server, extra config parameters.
#
# Calls the appropiate process: data processing pipeline or gem2s.
#
wrapper <- function(input_json) {
    # Get data from state machine input.
    input <- RJSONIO::fromJSON(input_json, simplify = FALSE)

    str(input)

    # common to gem2s and data processing
    task_name <- input$taskName
    server <- input$server
    input <- input[names(input) != "server"]
    pipeline_config <- load_config(server)
    process_name <- input$processName

    if (process_name == 'qc') {
        message_id <- call_data_processing(task_name, input, pipeline_config)
    } else if (process_name == 'gem2s') {
        message_id <- call_gem2s(task_name, input, pipeline_config)
    } else {
        stop("Process name not recognized.")
    }

    return(message_id)
}

#
# Init()
# Loads configuration and repeats the wrapper call until no more messages are received.
#
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
            withErrorTracing({
                message("Input ", input, " found")
                wrapper(input)

                message('Send task success')
                states$send_task_success(
                    taskToken = taskToken,
                    output = "{}"
                )
            }),
            error = function(e) {
                input_parse <- RJSONIO::fromJSON(input, simplify = FALSE)
                sample_text <- ifelse(is.null(input_parse$sampleUuid),
                                      "",
                                      paste0(" for sample ", input_parse$sampleUuid))

                error_txt <- paste0("R error at filter step ",
                                    input_parse$taskName, sample_text, "! : ", e$message)

                message(error_txt)
                states$send_task_failure(
                    taskToken = taskToken,
                    error = "We had an issue while processing your data.",
                    cause = error_txt
                )

                send_pipeline_fail_update(pipeline_config, input_parse$experimentId, input_parse$processName, "Error message placeholder")

                message("Sent task failure to state machine task: ", taskToken)
                message("recovered from error:", e$message)
            })
    }
}

init()
