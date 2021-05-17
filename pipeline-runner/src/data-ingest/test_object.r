library(testthat)
library(Seurat)

test_object <- function(scdata){ 

    test_that("Seurat object validation", {
        expect_is(scdata, "Seurat")
        expect_true(nrow(scdata) > 100, "Seurat") # more than 100 genes
        expect_true(scdata@active.assay == "RNA")
    })


    test_that("Validating metadata", {
        expect_is(scdata@meta.data, "data.frame")
        expect_true(nrow(scdata@meta.data) == ncol(scdata))
        md <- scdata@meta.data
        expect_true("barcode" %in% colnames(md))
        expect_true("orig.ident" %in% colnames(md))  
        if(any(grepl("^mt-", annotations$name, ignore.case = T))){
            expect_true("percent.mt" %in% colnames(md))
        }
        expect_true("doublet_scores" %in% colnames(md))
        expect_true("cells_id" %in% colnames(md))
        expect_true("samples" %in% colnames(md))

        test_that("Cell ids", {
            cellNumber <-ncol(scdata@assays$RNA@data)
            expect_equal(md$cells_id,0:(cellNumber-1))
        })  

        test_that("Percent mitocondrial", {
            expect_true(max(md$percent.mt) <= 100)
            expect_true(min(md$percent.mt) >= 0)
            #Polemic check to see that we are in percent and not fraction
            expect_true(max(md$percent.mt) > 1 || all(md$percent.mt==0))
        })  
    })


    test_that("Validating misc", {
        expect_is(scdata@misc, "list")
        misc <- scdata@misc
        expect_true("gene_annotations" %in% names(misc))
        expect_true("color_pool" %in% names(misc))
        expect_is(misc$color_pool,"character")
        expect_true(all(misc$gene_annotations$input == rownames(scdata)))
        expect_true(sum(duplicated(misc$gene_annotations$name))==0)
        #check that in duplicated positions (including the first) we have the gene id instead of the name.  
    })
}

