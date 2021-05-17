################################################
## STEP 4. Number of genes vs UMIs filter 
#################################################
# This filter focuses on filter cells that are far from the behaviour of the relationship between the number of genes (it measures the number of 
# genes in a cell that has at least one count) and the number of UMIs/molecules (the total number of counts in a cell). 
#'@Description Eliminates cells based on a p value and a linear regression generated from numGenes vs numUmis
#'@param Config
#    "numGenesVsNumUmis": {
#        "filterSettings": {
#            "regressionType": "gam",
#             "regressionTypeSettings": {
#                  "gam": {
#                       "p.level": 0.001
#                   }
#               }
#       },
#       "enabled": true
#       "auto": true
#   }
#

source('utils.r')

#' @description Filters seurat object based on classifier filter
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_numGenesVsNumUmis)
#'          - filterSettings: slot with thresholds
#'              - regressionType: String. Regression to be used: {gam}
#'              - regressionTypeSettings: list with the config settings for all the regression type options
#'                          - gam: for the gam option there is only one element: 
#'                                - p.level: which refers to  confidence level for deviation from the main trend
#' @export return a list with the filtered seurat object by numGenesVsNumUmis, the config and the plot values


task <- function(seurat_obj, config, task_name, sample_id, num_cells_to_downsample = 6000){
    options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
    print(paste("Running",task_name,sep=" "))
    print("Config:")
    print(config)
    print(paste0("Cells per sample before filter for sample ", sample_id))
    print(table(seurat_obj$orig.ident))
    #The format of the sample_id is
    # sample-WT1
    # we need to get only the last part, in order to grep the object.
    tmp_sample <- sub("sample-","",sample_id)

    # For now, we can get direcly p.level, but when we add more methods need to be change
    p.level <- config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]]$p.level

    # Check if it is required to compute sensible values. Sensible values are based on the funciton "gene.vs.molecule.cell.filter" from the pagoda2 package
    if (exists('auto', where=config)){
        if (as.logical(toupper(config$auto)))
            p.level <-  min(0.001, 1/ncol(seurat_obj))
    }

   # Check whether the filter is set to true or false
    if (as.logical(toupper(config$enabled))){
        # For now, we are going to suppor only gam as a linear model by robust estimation
        if (config$filterSettings$regressionType=="gam"){
            #Subsetting this sample
            obj_metadata <- seurat_obj@meta.data
            barcode_names_this_sample <- rownames(obj_metadata[grep(tmp_sample, rownames(obj_metadata)),]) 
            if(length(barcode_names_this_sample)==0){
                plots <- list()
                plots[generate_plotuuid(sample_id, task_name, 0)] <- list()
                return(list(data = seurat_obj,config = config,plotData = plots)) 
            }
            sample_subset <- subset(seurat_obj, cells = barcode_names_this_sample)
            # We regress the molecules vs the genes. This information are stored in nCount_RNA and nFeature_RNA respectively
            df <- data.frame(molecules = sample_subset$nCount_RNA, genes = sample_subset$nFeature_RNA)
            # We take log10 following the plot from the mock-up 
            df <- log10(df)
            # Rename the rows to be able to identify the valid cells
            rownames(df) <- colnames(sample_subset)
            df <- df[order(df$molecules, decreasing = FALSE), ]
            m <- MASS::rlm(genes ~ molecules, data = df)
            # Get the interval based on p.level paramter
            suppressWarnings(pb <- data.frame(predict(m, interval = "prediction", 
                level = 1 - p.level, type = "response")))
            # Define the outliers those that are below the lower confidence band and above the upper one.
            outliers <- rownames(df)[df$genes > pb$upr | df$genes < pb$lwr]

            # Keep the ones that are not outlier
            #Because we have subseted the object before creating the dataframes, and just selected the outlier values from that df
            #In the end we wil only remove values that are in the desired sample. 
            seurat_obj.filtered <- subset_safe(seurat_obj,colnames(seurat_obj)[!colnames(seurat_obj)%in%outliers])
                #Similarly, when creating the plot from the df, we will just be using the values for the required sample.
            plot1_data <- unname(purrr::map2(df$molecules,df$genes,function(x,y){c("log_molecules"=x,"log_genes"=y)}))
            plot1_data <- purrr::map2(plot1_data,unname(pb$lwr),function(x,y){append(x,c("lower_cutoff"=y))})
            plot1_data <- purrr::map2(plot1_data,unname(pb$upr),function(x,y){append(x,c("upper_cutoff"=y))})
            # have downsampling done
            
            # Downsample plotData
            # Handle when the number of remaining cells is less than the number of cells to downsample
            num_cells_to_downsample <- downsample_plotdata(ncol(sample_subset), num_cells_to_downsample)
            print(paste('sample of size', ncol(sample_subset), 'downsampled to', num_cells_to_downsample, 'cells'))
      
            set.seed(123)
            cells_position_to_keep <- sample(1:ncol(sample_subset), num_cells_to_downsample, replace = FALSE)
            cells_position_to_keep <- sort(cells_position_to_keep)
            plot1_data <- plot1_data[cells_position_to_keep]
        }
    }else{
        seurat_obj.filtered <- seurat_obj
        plot1_data <- list()
    }
    # update config
    config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]][["p.level"]] <- p.level



    # Scatter plot which is composed of:
    # x-axis: log_10_UMIs
    # y-axis: log_10_genes
    # bands that are conformed with the upper_cutoff and the lower_cutoff. We can print a band or dotted lines. 
    # Q: Should we return the point out the cells that are going to be excluded from the R side or this task can be done in 
    # the UI side.  
    plots <- list()
    plots[generate_plotuuid(sample_id, task_name, 0)] <- list(plot1_data)

    print(paste0("Cells per sample after filter for sample ", sample_id))
    print(table(seurat_obj.filtered$orig.ident))
    
    # the result object will have to conform to this format: {data, config, plotData : {plot1}}
    result <- list(
        data = seurat_obj.filtered,
        config = config,
        plotData = plots
    )

    return(result)
}


