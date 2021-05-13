require("RJSONIO")
require("paws")
require("zeallot")
require("ids")


source("handle_data.r")
source("utils.r")

# if not attached can cause errors when accessing metadata
require("Seurat")

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
    
    
    config[["aws_config"]] <- list(
        region = config$aws_region
    )
    
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
        assignInNamespace("update_endpoint_for_s3_config",function(request) {
            return(request)
        },ns="paws.common")
        
    }
    
    config[["originals_bucket"]] <- paste("biomage-originals", config$cluster_env, sep = "-")
    config[["source_bucket"]] <- paste("biomage-source", config$cluster_env, sep = "-")
    config[["processed_bucket"]] <- paste("processed-matrix", config$cluster_env, sep = "-")
    config[["results_bucket"]] <- paste("worker-results", config$cluster_env, sep = "-")
    config[["plot_data_bucket"]] <- paste("plots-tables", config$cluster_env, sep = "-")
    config[["sns_topic"]] <- paste("work-results", config$cluster_env, config$sandbox_id, sep = "-")
    config[["sns_topic"]] <- paste("arn:aws:sns", config$aws_region, config$aws_account_id, config$sns_topic, sep = ":")
    
    return(config)
}

#name changed from runstep
run_processing_step <- function(scdata, config, task_name, sample_id, debug_config) {
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
    handle_debug(scdata, config, task_name, sample_id, debug_config)
    out <- task(scdata, config, task_name, sample_id)
    return(out)
}

#
# input needs:
# sampleIds for all samples
# sampleNames for all samples
# projectId
# pipeline_config as defined by load_config
#
run_gem2s_step <- function(task_name,input,pipeline_config){
    switch(task_name,  
            download_gem = {
                import::here("/src/data-ingest/0-download_gem2s.r", task)
            },
            preproc = {
                import::here("/src/data-ingest/1_Preproc.r", task)
            },
            emptyDrops = {
                import::here("/src/data-ingest/2-1_Compute-metrics_emptyDrops.r",task)
            },
            doubletScores = {
                import::here("/src/data-ingest/2-2_Compute-metrics_doublets.r", task)
            },
            createSeurat = {
                import::here("/src/data-ingest/3_Seurat.r", task) 
            },
            prepareExperiment = {
                import::here("/src/data-ingest/4_Prepare_experiment.r", task)
            },
            stop(paste("Invalid task name given:", task_name))
    )
    task(input,pipeline_config)
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
    

    if (process_name == 'gem2s') {
        message_id <- run_gem2s_step(task_name,input,pipeline_config) 
    }else if(process_name == 'pipeline'){
        message_id <- call_data_processing(task_name, input, pipeline_config)
    }else{
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
