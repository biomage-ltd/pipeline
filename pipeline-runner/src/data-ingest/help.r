
# get_doublet_score function 
#' @description Get the cells with its doublet scores computed previously through scrublet
#' @param sample name of the sample to retrieve the doublet scores
#' 
#' @export data.frame with barcodes and doublet scores

get_doublet_score <- function(sample) {
    scores <-
        data.table::fread(
            paste("/output/doublet-scores-", sample, ".csv", sep = ""),
        )

    colnames(scores) <- c("barcodes", "doublet_scores", "doublet_class")
    rownames(scores) <- scores$barcodes    
    return(as.data.frame(scores))
}



# check_config function 
#' @description Create metadata dataframe from config files
#' @param scdata matrix with barcodes as columns to assign metadata information
#' @param sample name of the sample to retrieve the metadata
#' @param config config list from meta.json
#' 
#' @export dataframe with metadata information of the sample

check_config <- function(scdata, sample, config){
    metadata <- NULL
    metadata <- data.frame(row.names = colnames(scdata), samples=rep(sample, ncol(scdata)))

    # Check if "metadata" exists on config. If it is TRUE, we have other metadata information that we are
    # going to include in our experiment.
    if("metadata" %in% names(config)){
        rest_metadata <- as.data.frame(config$metadata)
        rest_metadata$sample <- ifelse(length(config$samples)>1, config$samples, sample)
        for(var in colnames(rest_metadata)){
            metadata[, var] <- rest_metadata[, var][match(metadata$sample, rest_metadata$sample)]
        }
    }
    
    return(metadata)
}