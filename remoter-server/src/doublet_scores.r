runDoubletScores <- function(data, config) {
    message("Filter 5: doubletScores")
    result <- doubletScores(data, config.doubletScores)
    return(result)
}

################################################
## STEP 5. Doublet score filter
#################################################
# This is a simplest filter that looks at a threshold value for the doublet scores.
# To separate cells with low droplet score from the ones that have a high droplet score content what makes us think that the are mistakenly considered as a single cell but they are actully two or more.
# This can be a useful first guess. The settings for such a filter can also contain a simple "probabilityThreshold" setting.


# The most uses values in doublet scores reporting in the scrublet paper [1] are around 0.25. There are not too much literature about how to compute
# a threshold. For now, we will offer two methods:
# --> Absolute threshold: In order to be not too extrictive the threshold is set to 0.25
generate_default_values_doubletScores <- function(scdata, config) {

    # HARDCODE
    threshold <- 0.25

    return(threshold)
}

#' @description Filters seurat object based on doubletScores scores
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_doubletScores)
#'          - filterSettings: slot with thresholds
#'                  - probabilityThreshold: Float. cut-off for the maximun probability scores that could have a cell. The doubletScores scores have been computed
#'                  through "scrublet" [1]
#'                  - binStep: Float. Bin size for the histogram
#' @export return a list with the filtered seurat object by doublet score, the config and the plot values

doubletScores <- function(scdata, config) {

    # Check if the experiment has doubletScores
    if (!"doublet_scores" %in% colnames(scdata@meta.data)) {
        message("Warning! No doubletScores scores has been computed for this experiment!")
        return(scdata)
    }

    probabilityThreshold <- config$filterSettings$probabilityThreshold

    # Check if it is required to compute sensible values. From the function 'generate_default_values_doubletScores', it is expected
    # to get a value --> probabilityThreshold.
    if (as.logical(toupper(config$auto))) {
        probabilityThreshold <- generate_default_values_doubletScores(scdata, config)
    }

    # Check wheter the filter is set to true or false
    if (as.logical(toupper(config$enabled))) {
        # Information regarding doublet score is pre-computed during the 'data-ingest'.
        scdata.filtered <- subset(scdata, subset = doublet_scores <= probabilityThreshold)
    } else {
        scdata.filtered <- scdata
    }

    # update config
    config$filterSettings$probabilityThreshold <- probabilityThreshold

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = scdata.filtered,
        config = config,
        plotData = list(
            # plot 1: histgram of doublet scores
            # AAACCCAAGCGCCCAT-1 AAACCCAAGGTTCCGC-1 AAACCCACAGAGTTGG-1
            #              0.161              0.198              0.284  ...
            doubletFilterHistogram = scdata$doublet_scores
        )
    )

    return(result)
}


# [1] Wolock SL, Lopez R, Klein AM. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data.
# Cell Syst. 2019 Apr 24;8(4):281-291.e9. doi: 10.1016/j.cels.2018.11.005. Epub 2019 Apr 3. PMID: 30954476; PMCID: PMC6625319.