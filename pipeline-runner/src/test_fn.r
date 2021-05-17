# This is a dummy function.
# Its only purpose is to demonstrate the input and output format
# that you will expect to your pipeline steps.

source('utils.r')

task <- function(input_data, input_config, task_name, sample_id) {
    # the result object will have to conform to this format.
    result <- list(
        data = input_data,
        config = input_config,
        plotData = list()
    )

    return(result)
}