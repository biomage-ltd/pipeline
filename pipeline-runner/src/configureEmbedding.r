################################################
## STEP 7. Compute embedding
#################################################
# Compute embedding step where we run dimensional reduction technniques such as t-SNE and  UMAP. Moreover, the cluster analysis is
# done also in this step. 

############ NEED TO CHECK CONFIG SCHEMA

# nPCs ?

#"configureEmbedding": {
#    "embeddingSettings": {
#        "method": "umap",
#        "methodSettings": {
#            "umap": {
#                "minimumDistance": 0.2,
#                "distanceMetric": "euclidean"
#            },
#            "tsne": {
#                "perplexity": 30,
#                "learningRate": 200
#            }
#        }
#    },
#    "clusteringSettings": {
#        "method": "louvain",
#        "methodSettings": {
#            "louvain": {
#                "resolution": 0.5
#            }
#        }
#    }
#},


source('utils.r')

task <- function(scdata, config, task_name, sample_id){

    # Check wheter the filter is set to true or false
    # Q: Can we disable the configureEmbedding?
    if (is.null(config$enabled) || as.logical(toupper(config$enabled))){
        scdata.embedding <- run_Embedding(scdata, config)
        scdata.embedding <- run_Clustering(scdata.embedding, config)
        scdata.embedding <- coloring_samples_and_cluster(scdata.embedding)

    }else{
        scdata.embedding <- scdata
    }

    cells_order <- rownames(scdata.embedding@meta.data)

    # Create a generic pattern with the UMAP coordinates
    embeddingPreview <- unname(purrr::map2(
        scdata.embedding@reductions$umap@cell.embeddings[cells_order, 1],
        scdata.embedding@reductions$umap@cell.embeddings[cells_order, 2],
        function(x,y) {c("x"=x,"y"=y)}
        )
    )

    # Create plot for cellsets
    embeddingPreviewByCellSets <- purrr::map2(embeddingPreview,
        unname(scdata.embedding@active.ident[cells_order]),
        function(x,y){append(x,c("cluster"=paste("Cluster", y, sep = " ")))}
    )

    # Adding color for cellsets
    embeddingPreviewByCellSets <- purrr::map2(embeddingPreviewByCellSets,
        unname(scdata.embedding@meta.data[cells_order, "color_active_ident"]),
        function(x,y){append(x,c("col"=y))}
    )

    # Create plot for samples
    embeddingPreviewBySamples <- purrr::map2(embeddingPreview,
        unname(scdata.embedding@meta.data[cells_order, "samples"]),
        function(x,y){append(x,c("sample"=y))}
    )
    
    # Adding color for samples
    embeddingPreviewBySamples <- purrr::map2(embeddingPreviewBySamples,
        unname(scdata.embedding@meta.data[cells_order, "color_samples"]),
        function(x,y){append(x,c("col"=y))}
    )

    # Create plot for MT-content
    embeddingPreviewMitochondrialContent <- purrr::map2(embeddingPreview,
        unname(scdata.embedding@meta.data[cells_order, "percent.mt"]),
        function(x,y){append(x,c("mt-content"=y))}
    )
    # Create plot for Doublet-score
    embeddingPreviewDoubletScore <- purrr::map2(embeddingPreview,
        unname(scdata.embedding@meta.data[cells_order, "doublet_scores"]),
        function(x,y){append(x,c("doublet-score"=y))}
    )

    plots <- list()
    plots[generate_plotuuid("", task_name, 0)] <- list(embeddingPreviewByCellSets)
    plots[generate_plotuuid("", task_name, 1)] <- list(embeddingPreviewBySamples)
    plots[generate_plotuuid("", task_name, 2)] <- list(embeddingPreviewMitochondrialContent)
    plots[generate_plotuuid("", task_name, 3)] <- list(embeddingPreviewDoubletScore)

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = scdata.embedding,
        config = config,
        plotData = plots
    )

    return(result)

}

# This function covers
#   - Compute embedding: t-SNE or UMAP
run_Embedding <- function(scdata, config){
    
    #################
    # Embedding part
    #################

    # The threshold was selected in the dataIntegration step
    # Until we update the rds in S3 where the numPCs is stored in misc, we need to handle it by HARDCODE to 30
    if ("numPCs" %in% names(scdata@misc))
        pca_nPCs <- scdata@misc[["numPCs"]]
    else
        pca_nPCs <- 30


    embed_method <- config$embeddingSettings$method
    embed_settings <- config$embeddingSettings$methodSettings[[embed_method]]

    if (embed_method == "umap"){
        message("Running embedding --> umap")
        # I get same embedding moving to Seurat default of 'uwot' with return.model = TRUE
        # I think it makes sense to do this as uwot-learn is being deprecated
        # seq_len avoids e.g. ndims=c(1,0) if user selects 0pcs (a bit more defensive)
        # I'm gonna add a ticket to make sure UI can't ask for 0 or 1 PCs as this will cause error anyways
        # please delete these comment
        scdata <- Seurat::RunUMAP(scdata, 
                                reduction = 'pca', 
                                dims = seq_len(pca_nPCs), 
                                verbose = FALSE, 
                                min.dist = embed_settings$minimumDistance, 
                                metric = embed_settings$distanceMetric, 
                                return.model = TRUE)
      
    } else if (embed_method == "tsne"){
      
        if (as.logical(toupper(config$auto))){
            perplexity <- min(30, ncol(scdata)/100)
            learningRate <- max(200, ncol(scdata)/12)
          
        }  else {
            perplexity <- embed_settings$perplexity
            learningRate <- embed_settings$learningRate
        }
      
        message("Running embedding --> tsne")
      
        scdata <- Seurat::RunTSNE(scdata,
                                reduction = 'pca',
                                dims = seq_len(pca_nPCs),
                                perplexity = perplexity,
                                learning.rate = learningRate)
      
    }

    return(scdata)

}

run_Clustering <- function(scdata, config){

    #################
    # Clustering part
    #################

    message("Running clustering")

    if(config$clusteringSettings$method=="louvain"){
        clustering_method <- 1 #"Louvain"
        clustering_resolution <- config$clusteringSettings$methodSettings[["louvain"]]["resolution"]

        # HARDCODE
        annoy.metric = "cosine"
        scdata <- Seurat::FindNeighbors(scdata, k.param = 20, annoy.metric = annoy.metric, verbose=FALSE)
        scdata <- Seurat::FindClusters(scdata, resolution=clustering_resolution, verbose = FALSE, algorithm = clustering_method)

    }

    return(scdata)
}


# This function return a seurat object with two new slot in meta.data
# - color_active_ident: since in the active ident should be the recent clusters we match a cluster to a color
# - color_samples: matching color and samples

coloring_samples_and_cluster <- function(scdata){

    ##########################
    # Coloring active ident
    ###########################
    if ("color_pool" %in% names(scdata@misc)){
        color_pool <- scdata@misc[['color_pool']]
    }else{ # THIS SHOULD BE REMOVE ONCE THE EXPERIMENT HAS BEEN UPDATED WITH THE NEW VERSION OF THE DATA-INGEST 
        color_pool <- c("#e377c2","#8c564b","#d62728","#2ca02c","#ff7f0e","#1f77b4","#f8e71c","#3957ff","#d3fe14",
        "#c9080a","#fec7f8","#0b7b3e","#0bf0e9","#c203c8","#fd9b39","#906407","#98ba7f","#fe6794","#10b0ff",
        "#ac7bff","#fee7c0","#964c63","#1da49c","#0ad811","#bbd9fd","#fe6cfe","#297192","#d1a09c","#78579e","#81ffad",
        "#739400","#ca6949","#d9bf01","#646a58","#d5097e","#bb73a9","#ccf6e9","#9cb4b6","#b6a7d4","#9e8c62","#6e83c8",
        "#01af64","#a71afd","#cfe589","#d4ccd1","#fd4109","#bf8f0e","#2f786e","#4ed1a5","#d8bb7d","#a54509","#6a9276",
        "#a4777a","#fc12c9","#606f15","#3cc4d9","#f31c4e","#73616f","#f097c6","#fc8772","#92a6fe","#875b44","#699ab3",
        "#94bc19","#7d5bf0","#d24dfe","#c85b74","#68ff57","#b62347","#994b91","#646b8c","#977ab4","#d694fd","#c4d5b5",
        "#fdc4bd","#1cae05","#7bd972","#e9700a","#d08f5d","#8bb9e1","#fde945","#a29d98","#1682fb","#9ad9e0","#d6cafe",
        "#8d8328","#b091a7","#647579","#1f8d11","#e7eafd","#b9660b","#a4a644","#fec24c","#b1168c","#188cc1","#7ab297",
        "#4468ae","#c949a6")

    }
    scdata$color_active_ident <- color_pool[as.numeric(scdata@active.ident)]

    ##########################
    # Coloring samples
    ###########################
    remaining.colors <- color_pool[-c(1:length(unique(scdata@meta.data$color_active_ident)))]
    if ("samples"%in%colnames(scdata@meta.data)) # In that case we are in multisample experiment
        scdata@meta.data[, "color_samples"] <- remaining.colors[as.numeric(as.factor(scdata$samples))]
    else
        scdata@meta.data[, "color_samples"] <- remaining.colors[1]

    return(scdata)
}


