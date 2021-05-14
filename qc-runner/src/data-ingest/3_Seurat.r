################################################
## 3_Seurat.r
##  - Create a seurat object per sample
##  - Adding emtpyDrops,  scrublet and MT-content
##  - Create a file with flag_filtered
################################################

suppressWarnings(library(Seurat))
suppressWarnings(library(Matrix))
suppressWarnings(library(dplyr))
suppressWarnings(require(data.table))
suppressWarnings(library(gprofiler2))
source("help.r")




################################################
## GETTING METADATA AND ANNOTATION
################################################

# adding_metrics_and_annotation function
#' @description We are going to process one sample in a seurat object. For it, we need to do:
#'  - Identify the possible metadata in the config file
#'  - Create a seurat object with the input of an sparse matrix
#'  - Annotate the gene in order to identify MT content in hsapiens and mmusculus
#'  - Computing MT-content
#'  - Getting scrublets
#'  - Getting emptyDrops
#'  - Identify flag_filtered
#'  - Save the rds with the seurat object for 1 sample
#' @param scdata Raw sparse matrix with the counts for one sample.
#' @param sample_name Name of the sample that we are preparing.
#' @param config Config of the project
#'
#' @return in the case that the input data was pre-filtered, we return a flag in order to disable the classifier filter. 
#' This flag is going to be store in the dynamoDB inside the samples-table.
adding_metrics_and_annotation <- function(scdata, sample, config, min.cells = 3, min.features = 10){
    message("Converting into seurat object sample --> ", sample)
    
    metadata <- check_config(scdata, sample, config)
    seurat_obj <- Seurat::CreateSeuratObject(scdata, assay='RNA', min.cells=min.cells, min.features=min.features, meta.data=metadata, project = config$name)
    
    organism <- config$organism

    # message("[", sample, "] \t finding genome annotations for genes...")
    #annotations <- gprofiler2::gconvert(
    #query = rownames(seurat_obj), organism = organism, target="ENSG", mthreshold = Inf, filter_na = FALSE)

    annotations <- read.delim("/output/features_annotations.tsv")

    if(any(grepl("^mt-", annotations$name, ignore.case = T))){
        message("[", sample, "] \t Adding MT information...")
        mt.features <-  annotations$input[grep("^mt-", annotations$name, ignore.case = T)]
        mt.features <- mt.features[mt.features %in% rownames(seurat_obj)]
        if (length(mt.features))
            seurat_obj <- PercentageFeatureSet(seurat_obj, features=mt.features , col.name = "percent.mt")
    }

    if (is.null(seurat_obj@meta.data$percent.mt)) seurat_obj$percent.mt <- 0

    message("[", sample, "] \t Getting scrublet results...")
    scores <- get_doublet_score(sample)
    rownames(scores) <- scores$barcodes

    idt <- scores$barcodes[scores$barcodes %in% rownames(seurat_obj@meta.data)]
    seurat_obj@meta.data[idt, "doublet_scores"] <- scores[idt, "doublet_scores"]
    # Doublet class is the classifications that scDblFinder does to set the threshold of doublet_scores
    # (https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/2_scDblFinder.html#thresholding-and-local-calibration)
    seurat_obj@meta.data[idt, "doublet_class"] <- scores[idt, "doublet_class"]

    message("[", sample, "] \t Adding emptyDrops...")
    file_ed <-  paste("/output/pre-emptydrops-", sample,".rds", sep = "")
    
    if (file.exists(file_ed)) {
        seurat_obj@tools$flag_filtered <- FALSE
        message("\t \t getting emptyDrops results...")
        emptydrops_out <- readRDS(file = file_ed)

        emptydrops_out_df <- emptydrops_out %>%
            as.data.frame() %>%
            rlang::set_names(~ paste0("emptyDrops_", .)) %>%
            tibble::rownames_to_column("barcode")
  
        # adding emptydrops data to meta.data for later filtering (using left join)
        meta.data <- seurat_obj@meta.data  %>%
            tibble::rownames_to_column("barcode") %>%
            dplyr::left_join(emptydrops_out_df)
            rownames(meta.data) <- meta.data$barcode  
  
        message("\t \t Adding emptyDrops scores information...")
        seurat_obj@meta.data <- meta.data
        # previously (before joining into meta.data), results were just dumped as a additional slot
        # leaving the code here in case bugs arise from above solution
        # seurat_obj@tools$CalculateEmptyDrops <- emptydrops_out
  
    } else {
        # Later on, when creating the config file, the enable will look the value of flag_filtered to   deactivate the classifier filter
        message("\t \t emptyDrops results not present, skipping...")
        seurat_obj@meta.data$emptyDrops_FDR <- NA
        seurat_obj@tools$flag_filtered <- TRUE
    }

    if (!dir.exists("/output/rds_samples")) 
        dir.create("/output/rds_samples")
    
    message("[", sample, "] Saving R object...")
    saveRDS(seurat_obj, file = paste("/output/rds_samples/", sample,".rds", sep = ""), compress = FALSE)

    return(seurat_obj@tools$flag_filtered)

}

################################################
## LOADING SPARSE MATRIX AND CONFIGURATION
################################################


task <- function(input,pipeline_config){
    # Loading the list with the raw sparse matrixs
    print("Starting")
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

    message("Creating Seurat Object...")
    flag_filtered <- sapply(names(scdata_list), function(sample_name) adding_metrics_and_annotation(scdata_list[[sample_name]], sample_name, config))

    # Since we need to store the flag_filtered in the dynamoDB, we need to persist in a file just to be read in 5_Upload-to-aws.py
    df_flag_filtered <- data.frame(samples=samples, flag_filtered=ifelse(flag_filtered, "Filtered", "Unfiltered"))
    write.table(df_flag_filtered, "/output/df_flag_filtered.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

    message("Step 3 completed.")
    print(list.files(paste("/output",sep = "/"),all.files=TRUE,full.names=TRUE,recursive=TRUE))
}

