################################################
## 4_Prepare_experiment.r
##  - Merging the samples for the current experiment
##  - Adding metadata: cellsId, color_pool and gene annotation
##  - Making QC plots
##  - Preparing dataProcess json file
################################################

suppressWarnings(library(Seurat))
suppressWarnings(library(Matrix))
suppressWarnings(library(dplyr))
suppressWarnings(require(data.table))
suppressWarnings(library(gprofiler2))

set.seed(123)
options(future.globals.maxSize= 1000 * 1024 ^ 2)
source("test_object.r")
source("help.r")


cellSizeDistribution_config <- function(seurat_obj, config) {
    import::here("cellSizeDistribution.r",generate_default_values_cellSizeDistribution)    
    minCellSize <- generate_default_values_cellSizeDistribution(seurat_obj,config,1e2)
    # update config
    config$filterSettings$minCellSize <- minCellSize

    return(config)
}

# There are some config parameters that depends on the data it-self. In this file we are going to create the functions
# that allow us to compute the best config parameter for Data Processing in the doubletScores step.

# To identify intelligently the treshold we are going to use the logic inside scDblFinder, which creates a classification
# (singlet our doublet) [ref: https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/2_scDblFinder.html#thresholding-and-local-calibration]
# To set the auto value we are going to use as a threshold the maximun score that is given to a singlet. 

doubletScores_config <- function(scdata, config){

    # Minimun score that has a singlet 
    probabilityThreshold <-  max(scdata$doublet_scores[scdata$doublet_class=="singlet"])
    # update config
    config$filterSettings$probabilityThreshold <- probabilityThreshold

    return(config)
}



# There are some config parameters that depends on the data it-self. In this file we are going to create the functions
# that allow us to compute the best config parameter for Data Processing in the numGenesVsNumUmis step.

numGenesVsNumUmis_config <- function(scdata, config){

    # Sensible values are based on the funciton "gene.vs.molecule.cell.filter" from the pagoda2 package
    p.level <-  min(0.001, 1/ncol(scdata))
    # update config
    config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]]$p.level <- p.level

    return(config)
}



################################################
## SAVING CONFIG FILE 
################################################
# 
# We are going to store the final config to config_dataProcessing.json, in order to upload to dynamoDB.
# The unisample experiments does not require any change, but for the multisample experiment we need
# to add the filtering parameter for each sample (only in the steps that is required.)
# We are going to differentiate in samples only in the steps:
# --> cellSizeDistribution
# --> numGenesVsNumUmis
# --> doubletScores
#
# For both of them, we will run again the step fn for each sample (samples names are stored in metadata type)

# Function to recompute the step fn and store the new config of each sample inside the latest config file 
# We need to iterate per sample and compute separately the step fn.
# Example of structure:
# {
# “filterSettings”: {
#   “probabilityThreshold”: 0.2,
#   “binStep”: 0.05
# },
# “sample-KO”: {
#   “filterSettings”: {
#     “probabilityThreshold”: 0.1,
#     “binStep”: 100
#   }
# },
# “sample-WT1": {
#                     “filterSettings”: {
#                         “probabilityThreshold”: 0.1,
#                         “binStep”: 45
#                     }
#                 }
# }

add_custom_config_per_sample <- function(step_fn, config, scdata, samples){
  
  # We upadte the config file, so to be able to access the raw config we create a copy
  config.raw <- config
  
  for(sample in samples){
    # Downsample the seurat object to a unisample experiment
    scdata_sample <- subset(scdata, samples %in% sample)
    # Run the step fun with the unisample experiment and keep the config result
    result_config <- step_fn(scdata_sample, config.raw)
    # Inside the config of the samples we are not storing the auto and enable settings, so we remove them
    result_config$auto <- NULL
    result_config$enabled <- NULL
    # Update config with the unisample thresholds
    config[[paste("sample-", sample, sep = "")]] <- result_config
  }
  
  return(config)
  
}


task <- function(input,pipeline_config){
    ################################################
    ## Merging Seurat object
    ################################################

    message("Loading configuration...")
    config <- RJSONIO::fromJSON("/input/meta.json")

    # Check which samples have been selected. Otherwiser we are going to use all of them. 
    if (length(config$samples)>0){
        samples <- config$samples
    }else{
        samples <- gsub("\\..*", "", list.files("/output/rds_samples"))
    }


    message("Reloading samples rds for current experiment...")
    seurat_obj_list <- list()
    for (sample in samples){
        seurat_obj_list[[sample]] <- readRDS(paste("/output/rds_samples/", sample,".rds", sep = ""))
    }

    # Merging samples and adding a prefix with the sample name. In pipeline we grep in barcodes to filter by sample.
    if (length(seurat_obj_list)==1){
        seurat_obj <- seurat_obj_list[[1]]
        seurat_obj <- RenameCells(object = seurat_obj, add.cell.id = names(seurat_obj_list)[1])
    }else{
        seurat_obj <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1], add.cell.ids = c(samples))
    }


    ################################################
    ## Adding metadata
    ################################################

    message("Storing gene annotations...")
    organism <- config$organism
    #annotations <- gprofiler2::gconvert(
    #query = rownames(seurat_obj), organism = organism, target="ENSG", mthreshold = Inf, filter_na = FALSE)
    annotations <- read.delim("/output/features_annotations.tsv")

    # In order to avoid duplicated genes names, we are going to add the ENSEMBL ID for those 
    # genes that are duplicated (geneNameDuplicated-ENSEMBL)
    gname <- annotations$name
    # Keep original name in 'original_name' variable
    annotations$original_name <- gname
    is.dup <- duplicated(gname) | duplicated(gname, fromLast=TRUE)
    annotations$name[is.dup] <- paste(gname[is.dup], annotations$input[is.dup], sep = " - ")

    # Ensure index by rownames in seurat_obj
    annotations <- annotations[match(rownames(seurat_obj), annotations$input), ]
    rownames(annotations) <- annotations$input

    seurat_obj@misc[["gene_annotations"]] <- annotations

    message("Storing cells id...")
    # Keeping old version of ids starting from 0
    seurat_obj$cells_id <- 0:(nrow(seurat_obj@meta.data)-1)

    message("Storing color pool...")
    # We store the color pool in a slot in order to be able to access it during configureEmbedding
    color_pool <- RJSONIO::fromJSON("data-ingest/color_pool.json")
    seurat_obj@misc[["color_pool"]] <- color_pool
    message("Stored pool")
    
    ################################################
    ## Checking filtered data
    ################################################

    df_flag_filtered <- read.delim("/output/df_flag_filtered.txt")
    any_filtered <- "Filtered" %in% df_flag_filtered$flag_filtered
    message("saved filtered flag")
    # Im leaving the qc plots for now. If we dont use them we might want to remove them.
    ################################################
    ## QC plots
    ################################################

    # if (!dir.create("/output/QC_plots")) 
    #     dir.create("/output/QC_plots")

    # pdf("/output/QC_plots/mitochondrialFractionLogHistogram.pdf")
    # hist(seurat_obj$percent.mt, breaks=200)
    # dev.off()

    # pdf("/output/QC_plots/mitochondrialFractionLogScatter.pdf")
    # plot(seurat_obj$nCount_RNA, seurat_obj$percent.mt)
    # dev.off()

    # pdf("/output/QC_plots/UMI_hist.pdf",width = 14, height = 7)
    # par(mfcol=c(1,2),cex.lab=2)
    # hist(colSums(seurat_obj@assays$RNA@data), breaks=200)
    # hist(rowSums(seurat_obj@assays$RNA@data), breaks=200)
    # dev.off()

    # genes_umi <- rowSums(seurat_obj@assays$RNA@data)
    # counts <- colSums(seurat_obj@assays$RNA@data)

    # pdf("/output/QC_plots/nCount_vs_mito.fr_hist.pdf")
    # plot(seurat_obj@meta.data$nCount_RNA, seurat_obj@meta.data$percent.mt)
    # dev.off()

    # pdf("/output/QC_plots/GenesVsNumUmis.pdf")
    # plot(seurat_obj@meta.data$nCount_RNA, seurat_obj@meta.data$nFeature_RNA)
    # dev.off()


    # sink("/output/QC_plots/routput.Rout")
    # print("UMI counts of highest expressed genes")
    # tail(sort(genes_umi), n=30)
    # print("median rowsums")
    # median(rowSums(seurat_obj@assays$RNA@data))
    # sink()

    # # some dignostic plots for empty-drops
    # meta.data <- seurat_obj@meta.data
    # if(all(!is.na(meta.data$emptyDrops_FDR))){
    #     pdf("/output/QC_plots/emptyDrops_Total_hist.pdf",width = 16, height = 14)
    #     par(mfcol=c(2,2),cex.lab=2)
    #     plot1_data_1  <- log10(meta.data$emptyDrops_FDR)
    #     names(plot1_data_1) <- rep("FDR", length(plot1_data_1))
    #     plot1_data_2  <- log10(meta.data$nCount_RNA)
    #     names(plot1_data_2) <- rep("log_u", length(plot1_data_2))
    #     plot(log10(meta.data$nCount_RNA),log10(meta.data$emptyDrops_FDR))
    #     plot((meta.data$nCount_RNA),(meta.data$emptyDrops_FDR))
    #     dev.off()
    # }

    ################################################
    ## Testing Seurat object before save
    ################################################

    #test_object(seurat_obj)

    ################################################
    ## Saving files
    ################################################


    message("saving R object...")
    saveRDS(seurat_obj, file = "/output/experiment.rds", compress = FALSE)





    message("saving multsiample info...")
    write.table(
        data.frame(Cells_ID=seurat_obj$cells_id, Value=seurat_obj$samples),
        file = "/output/samples-cells.csv",
        quote = F, col.names = F, row.names = F,
        sep = "\t"
    )


    if("metadata" %in% names(config)){
        variables_metadata <- names(config$metadata)
        metadata_dynamo <- seurat_obj@meta.data[, c("cells_id", variables_metadata)]

        message("saving multsiample info...")
        write.table(
            metadata_dynamo,
            file = "/output/metadata-cells.csv",
            quote = F, col.names = T, row.names = F,
            sep = "\t"
        )
    }

    write.table(
        colnames(seurat_obj),
        file = "/output/r-out-cells.csv",
        quote = F, col.names = F, row.names = F,
        sep = "\t"
    )

    write.table(
        seurat_obj@misc[["gene_annotations"]][seurat_obj@misc[["gene_annotations"]]$input%in%rownames(seurat_obj), ],
        file = "/output/r-out-annotations.csv",
        quote = F, col.names = F, row.names = F,
        sep = "\t"
    )

    message("saving normalized matrix...")
    #Matrix::writeMM(t(seurat_obj@assays[[seurat_obj@active.assay]]@data
    #                ), file = "/output/r-out-normalized.mtx")

    #message("saving raw matrix...")
    #Matrix::writeMM(t(
    #    seurat_obj@assays[["RNA"]]@counts[
    #        rownames(seurat_obj),
    #        colnames(seurat_obj)
    #    ]),
    #    file = "/output/r-out-raw.mtx", verb
    #)
    print(seurat_obj)


    ################################################
    ## DATA PROCESSING
    ################################################

    #[HARDCODED]
    config.cellSizeDistribution <- list(enabled="true", 
        auto="true", 
        filterSettings = list(minCellSize=1080, binStep = 200)
    )

    config.mitochondrialContent <- list(enabled="true", auto="true", 
        filterSettings = list(method="absolute_threshold", methodSettings = list(
            absolute_threshold=list(maxFraction=0.1, binStep=0.05)
            )
        )
    )

    config.classifier <- list(enabled=tolower(as.character(!any_filtered)) # emptyDrops results not present
        , auto="true", 
        filterSettings = list(FDR=0.01)
    )

    config.numGenesVsNumUmis <- list(enabled="true", auto="true", 
        filterSettings = list(regressionType = "gam", regressionTypeSettings = list(
            "gam" = list(p.level=0.001)
            )
        )
    )

    config.doubletScores <- list(enabled="true", auto="true", 
        filterSettings = list(probabilityThreshold = 0.5, binStep = 0.05)
    )

    # BE CAREFUL! The method is based on config.json. For multisample only seuratv4, for unisample LogNormalize
    # hardcoded because unisample check is performed in dataIntegration 
    identified.method <- 'seuratv4'
    config.dataIntegration <- list(auto="false", 
        dataIntegration = list( method = identified.method , 
                            methodSettings = list(seuratv4=list(numGenes=2000, normalisation="logNormalize"), 
                                                unisample=list(numGenes=2000, normalisation="logNormalize"))),
        dimensionalityReduction = list(method = "rpca", numPCs = 30, excludeGeneCategories = c())
    )

    config.configureEmbedding <- list(auto="false", 
        embeddingSettings = list(method = "umap", methodSettings = list(
                                    umap = list(minimumDistance=0.3, distanceMetric="euclidean"), 
                                    tsne = list(perplexity=min(30, ncol(seurat_obj)/100), learningRate=max(200, ncol(seurat_obj)/12))
                                ) 
                            ), 
        clusteringSettings = list(method = "louvain", methodSettings = list(
                                louvain = list(resolution = 0.5)
                                )
                            )
    )

    # Compute for multisample and unisample
    config.cellSizeDistribution <- add_custom_config_per_sample(cellSizeDistribution_config, config.cellSizeDistribution, seurat_obj, unique(seurat_obj$samples))
    config.numGenesVsNumUmis <- add_custom_config_per_sample(numGenesVsNumUmis_config, config.numGenesVsNumUmis, seurat_obj, unique(seurat_obj$samples))
    config.doubletScores <- add_custom_config_per_sample(doubletScores_config, config.doubletScores, seurat_obj, unique(seurat_obj$samples))

    # When we remove the steps from data-ingest we need to change here the default config. 
    # Save config for all steps. 
    config <- list(
    cellSizeDistribution = config.cellSizeDistribution
    , mitochondrialContent = config.mitochondrialContent
    , classifier = config.classifier
    , numGenesVsNumUmis = config.numGenesVsNumUmis
    , doubletScores = config.doubletScores
    , dataIntegration = config.dataIntegration
    , configureEmbedding = config.configureEmbedding
    )


    # Export to json
    exportJson <- RJSONIO::toJSON(config, pretty = T)
    # The RJSONIO library add '' to boolean keys, so we will remove them.
    exportJson <- gsub('\"true\"', "true", exportJson)
    exportJson <- gsub('\"false\"', "false", exportJson)
    # Trnasform null into []
    exportJson <- gsub('null', "[]", exportJson)
    message("config file...")
    write(exportJson, "/output/config_dataProcessing.json")

    message("Step 4 completed.")
    print(list.files(paste("/output",sep = "/"),all.files=TRUE,full.names=TRUE,recursive=TRUE))
}


