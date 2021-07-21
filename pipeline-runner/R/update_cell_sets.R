update_cell_sets <- function(scdata, experiment_id, pipeline_config) {

    # Design cell_set scratchpad for DynamoDB
    scratchpad <- list(
        key = "scratchpad",
        name = "Scratchpad",
        rootNode = TRUE,
        children = list(),
        type = "cellSets"
    )

    color_pool <- get_color_pool()
    samples_set <- samples_sets(scdata, color_pool)

    # remove used colors from pool
    used <- seq_along(samples_set$children)
    color_pool <- color_pool[-used]

    # Design cell_set meta_data for DynamoDB
    cell_sets <- c(list(scratchpad), list(samples_set))

    if ("meta_vars" %in% names(scdata@misc)) {
        cell_sets <- c(cell_sets, meta_sets(scdata, color_pool))
    }

    cell_sets <- list(cellSets = cell_sets)

    # cell sets file to s3
    cell_sets_data <- RJSONIO::toJSON(cell_sets)

    put_object_in_s3(pipeline_config,
                     bucket = pipeline_config$cell_sets_bucket,
                     object = charToRaw(cell_sets_data),
                     key = experiment_id
    )

}

samples_sets <- function(scdata, color_pool) {

    annot <- data.frame(cell_id = scdata$cells_id, sample_id = scdata$samples)

    cell_set <- list(
        key = "sample",
        name = "Samples",
        rootNode = TRUE,
        children = list(),
        type = "metadataCategorical"
    )

    sample_ids <- unique(scdata$samples)
    sample_names <- unique(scdata$sample_name)

    for (i in seq_along(sample_ids)) {
        sample_id <- toString(sample_ids[i])
        cell_ids <- annot$cell_id[annot$sample_id == sample_id]

        cell_set$children[[i]] <- list(
            key = sample_id,
            name = sample_names[i],
            color = color_pool[i],
            cellIds = cell_ids
        )
    }

    return(cell_set)
}


# cell_sets fn for seurat metadata information
meta_sets <- function(scdata, color_pool) {

    meta_vars <- scdata@misc[['meta_vars']]
    cell_set_list <- c()

    # The first column is the cells_id, the rest is the metadata information
    for (meta_var in meta_vars) {

        cell_set <- list(
            "key" = meta_var,
            "name" = meta_var,
            "rootNode" = TRUE,
            "children" = c(),
            "type" = "metadataCategorical"
        )

        annot <- scdata@meta.data[[meta_var]]
        values <- unique(annot)

        for (i in seq_along(values)) {
            value <- values[i]
            is.val <- annot == value
            cell_ids <- unname(scdata$cells_id[is.val])


            cell_set$children[[i]] <- list(
                "key" = paste(meta_var, value, sep = "-"),
                "name" = value,
                "color" = color_pool[i],
                "cellIds" = cell_ids
            )
        }
        cell_set_list <- c(cell_set_list, list(cell_set))
    }
    return(cell_set_list)
}

