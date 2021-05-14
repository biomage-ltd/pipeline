################################################
## 2-2_Compute-metrics_scrublet.r
##  - Compute doublet scores per sample
##  - Save an txt file with doublet scores per sample [doublet-scores-sampleName]
################################################

suppressWarnings(library(Seurat))
suppressWarnings(library(scDblFinder))


# compute_doublet_scores function
#' @description Save the result of doublets scores per sample. 
#' @param scdata Raw sparse matrix with the counts for one sample.
#' @param sample_name Name of the sample that we are preparing.
#'
#' @return 
compute_doublet_scores <- function(scdata, sample_name, min.features = 10) {
    message("loading scd")
    library(scDblFinder)
    message("loaded")
    message("Sample --> ", sample_name, "...")
    scdata_DS <- scDblFinder(scdata[, Matrix::colSums(scdata>0)>=min.features], dbr = NULL, trajectoryMode = FALSE)
    df_doublet_scores <- data.frame(Barcodes=rownames(scdata_DS@colData), doublet_scores=scdata_DS@colData$scDblFinder.score,
    doublet_class = scdata_DS@colData$scDblFinder.class)

    write.table(df_doublet_scores, file = paste("/output/doublet-scores-", sample_name, ".csv", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}



task <- function(input,pipeline_config){
    print(list.files(paste("/output",sep = "/"),all.files=TRUE,full.names=TRUE,recursive=TRUE))
    message("reloading old matrices...")
    
    scdata_list <- readRDS("/output/pre-doublet-scdata_list.rds")

    message("Loading configuration...")
    config <- RJSONIO::fromJSON("/input/meta.json")

    # Check which samples have been selected. Otherwiser we are going to use all of them. 
    if (length(config$samples)>0){
        samples <- config$samples
    }else{
        samples <- names(scdata_list)
    }

    scdata_list <- scdata_list[samples]

    message("calculating probability of barcodes being doublets...")
    for (sample_name in names(scdata_list)) {
        compute_doublet_scores(scdata_list[[sample_name]], sample_name)
    }

    message("Step 2-2 completed.")
    print(list.files(paste("/output",sep = "/"),all.files=TRUE,full.names=TRUE,recursive=TRUE))
}

