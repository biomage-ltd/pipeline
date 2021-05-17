# Biomage Pipeline

The Cellscope pipeline project for dependency-managed work processing.

## Getting started

The steps of the pipeline that are run through this project are started
spontaneously on your machine as Docker containers, simulating Kubernetes
in local development.

We have included a utility so you can automatically monitor containers spawned
and read their logs as they are executing.

For local development, you should already have Docker and Node.js installed, as well as
Inframock running.

Afterwards, you can install the pipeline dependencies with:

```bash
make install
```

To build and run the pipeline:

```bash
make build && make run
```

A similar message should appear:

```
> node src/app.js

Loading CloudFormation for local container launcher...
Creating mock Lambda function on InfraMock...
No previous stack found on InfraMock.
Stack with ARN arn:aws:cloudformation:eu-west-1:000000000000:stack/local-container-launcher/106d1df9 successfully created.
Waiting for Docker events...
```

Logs from pipelines run through the API will apear here.

## Rebuiling the docker images

```bash
make build
```

## Debugging locally
To save the parameters (`config`, `seurat_obj`, etc) to a task function, specify DEBUG_STEP and DEBUG_PATH.
Available tasks include all those defined in `run_step` of  [init.r](pipeline-runner/src/init.r) as well as `DEBUG_STEP=all` 
to save the parameters to all task functions:

```bash
# e.g. DEBUG_STEP=dataIntegration
DEBUG_STEP=task_name DEBUG_PATH=${PWD}/mydebug make run
```
Note: `DEBUG_PATH` needs to be an absolute path. The above can be used to populate the local subfolder `mydebug`.

When the pipeline is run, it will save the parameters to the specified `task_name` in `DEBUG_PATH`. You
can load these into your R environment:

```R
# clicking the file in RStudio does this for you
load('/path/to/debug/folder/{task_name}_{sample_id}.RData')

# if you need to load multiple tasks, you can load each into a seperate environment
# you would when access objects using e.g. task_env$scdata
task_env <- new.env()
load('/path/to/debug/folder/{task_name}_{sample_id}.RData', envir = task_env)
```
