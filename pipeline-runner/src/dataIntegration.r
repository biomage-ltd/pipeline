################################################
## STEP 6. Data integration
#################################################
# Data integration step where batch effect is corrected through data integration methods such as "Seurat V3". 
# The data integration step include the normalization and the PCA analysis.

# NEED TO DISCUSS: functio to compute default parameters

############ NEED TO CHECK CONFIG SCHEMA

#   "dataIntegration": {
#       "dataIntegration": {
#           "method": "seuratv4",
#           "methodSettings": {
#               "seuratv4": {
#                   "numGenes": 2000,
#                   "normalization": "logNormalise"
#               }
#           }
#       },
#       "dimensionalityReduction": {
#           "method": "rpca",  --> NOT USE!
#           "numPCs": 30,
#           "excludeGeneCategories": []
#       }
#   },

task <- function(scdata, config,task_name,sample_id){
    # this options shows the callstack when an error is thrown
    options(error = function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
    # increase maxSize from the default of 500MB to 32GB
    # TODO: ask Marcell for his opinion
    options(future.globals.maxSize= 32 * 1024 * 1024^2)
    # Check wheter the filter is set to true or false
    # So far we only support Seurat V3
    scdata.integrated <- run_dataIntegration(scdata, config)
    # Compute explained variance for the plot2
    eigValues = (scdata.integrated@reductions$pca@stdev)^2  ## EigenValues
    varExplained = eigValues / sum(eigValues)
    # As a short solution, we are going to store an intermediate slot for the numPCs, since this parameter is required when performing
    # the computeEmdedding. The main reason to do not have in the config.configureEmbedding is that this parameter does not change in the configureEmbedding step.
    scdata.integrated@misc[["numPCs"]] <- config$dimensionalityReduction$numPCs

    scdata.integrated <- colorObject(scdata.integrated)
    cells_order <- rownames(scdata.integrated@meta.data)
    plot1_data <- unname(purrr::map2(scdata.integrated@reductions$umap@cell.embeddings[, 1],scdata.integrated@reductions$umap@cell.embeddings[, 2],function(x,y){c("x"=x,"y"=y)}))
    
    #Adding color and sample id
    plot1_data <- purrr::map2(plot1_data,
        unname(scdata.integrated@meta.data[cells_order, "samples"]),
        function(x,y){append(x,list("sample"=y))}
    )
    plot1_data <- purrr::map2(plot1_data,
        unname(scdata.integrated@meta.data[cells_order, "color_samples"]),
        function(x,y){append(x,list("col"=y))}
    )

    plot2_data <- unname(purrr::map2(1:50,varExplained,function(x,y){c("PC"=x,"percentVariance"=y)}))

    plots <- list()
    plots[generate_plotuuid("", task_name, 0)] <- list(plot1_data)
    plots[generate_plotuuid("", task_name, 1)] <- list(plot2_data)    

    # For now config is not updated, since there is not new changes
    # config <- ...

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = scdata.integrated,
        config = config,
        plotData = plots
    )

    return(result)

}

# This function covers
#   - Integrate the data using the variable "type" (in case of data integration method is selected) and normalize using LogNormalize method.
#   - Compute PCA analysis
#   - To visualize the results of the batch effect, an UMAP with default setting has been made. 
# Seurat V3 pipeline (see for other methods: https://satijalab.org/seurat/archive/v3.0/integration.html)
run_dataIntegration <- function(scdata, config){
    
    method <- config$dataIntegration$method
    nfeatures <- config$dataIntegration$methodSettings[[method]]$numGenes
    normalization <- config$dataIntegration$methodSettings[[method]]$normalisation
    
    #Caps lock independent just in case anything goes wrong in dynamo or wherever. 
    if (toupper(normalization)=="LOGNORMALIZE"  | toupper(normalization)=="LOGNORMALISE") normalization<-"LogNormalize"
    # Q: Should be the numPCs input for the RunPCA? Since we are looking into the explained variance I suggest to RunPCA with 50 and only
    # use this parameter for the data integration purpose.
    numPCs <- config$dimensionalityReduction$numPCs

    #HARDCODE
    # We require just to get an overview of the data-integration
    umap_min_distance <- 0.3
    umap_distance_metric <- "euclidean"
    # default for FindIntegrationAnchors
    k.filter <- 200 

    # temporary to make sure we don't run integration if unisample
    nsamples <- length(unique(scdata$samples))
    
    # Currently, we only support Seurat V4 pipeline for the multisample integration
    if(nsamples > 1 && method=="seuratv4"){
        #FIX FOR CURRENT DATASET!!!!!!
        Seurat::DefaultAssay(scdata) <- "RNA"
        data.split <- Seurat::SplitObject(scdata, split.by = "samples")
        for (i in 1:length(data.split)) {
            data.split[[i]] <- Seurat::NormalizeData(data.split[[i]], normalization.method = normalization, verbose = F)
            data.split[[i]] <- Seurat::FindVariableFeatures(data.split[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
        }
        # If Number of anchor cells is less than k.filter/2, there is likely to be an error:
        # Note that this is a heuristic and was found to still fail for small data-sets

        # Try to integrate data (catch error most likely caused by too few cells)

        tryCatch({
          k.filter <- min(ceiling(sapply(data.split, ncol)/2), k.filter)
          data.anchors <- Seurat::FindIntegrationAnchors(object.list = data.split, dims = 1:numPCs, k.filter = k.filter, verbose = TRUE)

          # @misc slots not preserved so transfer
          misc <- scdata@misc
          scdata <- Seurat::IntegrateData(anchorset = data.anchors, dims = 1:numPCs)
          scdata@misc <- misc
          Seurat::DefaultAssay(scdata) <- "integrated"
        }, error = function(e){          # Specifying error message
          # ideally this should be passed to the UI as a error message:
          print(table(scdata$samples))
          print(e)
          print(paste("current k.filter:", k.filter))
          # Should we still continue if data is not integrated? No, right now..
          print("Current number of cells per sample: ")
          print(table(scdata$samples))
          warning("Error thrown in IntegrateData: Probably one/many of the samples contain to few cells.\nRule of thumb is that this can happen at around < 100 cells.")
          # An ideal solution would be to launch an error to the UI, howerver, for now, we will skip the integration method. 
          print("Skipping integration step")
          scdata <- Seurat::NormalizeData(scdata, normalization.method = normalization, verbose = F)

        })

    }else{
        print('Only one sample detected.')
        # Else, we are in unisample experiment and we only need to normalize 
        scdata <- Seurat::NormalizeData(scdata, normalization.method = normalization, verbose = F)
    }

    scdata <- FindVariableFeatures(scdata, selection.method = "vst", assay = "RNA", nfeatures = nfeatures, verbose = FALSE)
    vars <- HVFInfo(object = scdata, assay = "RNA", selection.method = 'vst') # to create vars
    annotations <- scdata@misc[["gene_annotations"]]
    vars$SYMBOL <- annotations$name[match(rownames(vars), annotations$input)]
    vars$ENSEMBL <- rownames(vars)
    scdata@misc[["gene_dispersion"]] <- vars

    # Scale in order to compute PCA
    scdata <- Seurat::ScaleData(scdata, verbose = F)

    # HARDCODE numPCs to 50
    scdata <- Seurat::RunPCA(scdata, npcs = 50, features = Seurat::VariableFeatures(object=scdata), verbose=FALSE)

    # Compute embedding with default setting to get an overview of the performance of the bath correction
    scdata <- Seurat::RunUMAP(scdata, reduction='pca', dims = 1:numPCs, verbose = F, umap.method = "uwot-learn", min.dist = umap_min_distance, metric = umap_distance_metric)

    return(scdata)
}

colorObject <- function(data){

    if ("color_pool" %in% names(data@misc)){
        color_pool <- data@misc[['color_pool']]
    }else{ # THIS SHOULD BE REMOVE ONCE THE EXPERIMENT HAS BEEN UPDATED WITH THE NEW VERSION OF THE DATA-INGEST 
        color_pool <- c("#e377c2","#8c564b","#d62728","#2ca02c","#ff7f0e","#1f77b4","#f8e71c","#3957ff","#d3fe14",
        "#c9080a","#fec7f8","#0b7b3e","#0bf0e9","#c203c8","#fd9b39","#888593","#906407","#98ba7f","#fe6794","#10b0ff",
        "#ac7bff","#fee7c0","#964c63","#1da49c","#0ad811","#bbd9fd","#fe6cfe","#297192","#d1a09c","#78579e","#81ffad",
        "#739400","#ca6949","#d9bf01","#646a58","#d5097e","#bb73a9","#ccf6e9","#9cb4b6","#b6a7d4","#9e8c62","#6e83c8",
        "#01af64","#a71afd","#cfe589","#d4ccd1","#fd4109","#bf8f0e","#2f786e","#4ed1a5","#d8bb7d","#a54509","#6a9276",
        "#a4777a","#fc12c9","#606f15","#3cc4d9","#f31c4e","#73616f","#f097c6","#fc8772","#92a6fe","#875b44","#699ab3",
        "#94bc19","#7d5bf0","#d24dfe","#c85b74","#68ff57","#b62347","#994b91","#646b8c","#977ab4","#d694fd","#c4d5b5",
        "#fdc4bd","#1cae05","#7bd972","#e9700a","#d08f5d","#8bb9e1","#fde945","#a29d98","#1682fb","#9ad9e0","#d6cafe",
        "#8d8328","#b091a7","#647579","#1f8d11","#e7eafd","#b9660b","#a4a644","#fec24c","#b1168c","#188cc1","#7ab297",
        "#4468ae","#c949a6")

    }
    data$color_active_ident <- color_pool[as.numeric(data@active.ident)]

    ##########################
    # Coloring samples
    ###########################
    if ("samples"%in%colnames(scdata@meta.data)) # In that case we are in multisample experiment
        data@meta.data[, "color_samples"] <- color_pool[as.numeric(as.factor(data$samples))]
    else
        data@meta.data[, "color_samples"] <- color_pool[1]
    return(data)
}
