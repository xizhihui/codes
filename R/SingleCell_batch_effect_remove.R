# prepare packages ----
install_pkgs <- function() {
    # pkg - batch correction ----
    devtools::install_github("immunogenomics/harmony")
    devtools::install_github('MacoskoLab/liger')
    BiocManager::install("SeuratWrappers")
    BiocManager::install("scran") # fastMNN
    BiocManager::install("scMerge")# need libgsl0-dev
    BiocManager::install("zinbwave")# need libgsl0-dev too
    install.packages('Seurat')
    library(reticulate)
    #use_python("/usr/bin/python")#use your own python
    py_install(c("scanorama", "bbknn"))
    
    # pkg - visualization ----
    BiocManager::install("clustree")
    BiocManager::install("riverplot")
    BiocManager::install("ggalluvial")
}
#install_pkgs()

# wrappers ----
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(liger)
    library(harmony)
    library(scran)
    library(scater)
    library(scMerge)
    library(zinbwave)
    library(Seurat)
    library(SeuratWrappers)
})
#reticulate::import(c("scanorama", "bbknn"))

`%>%` <- magrittr::`%>%`
split_SCE <- function(sce, split.by) {
    # return list of sce with assays and coldata,rowdata
    stopifnot(split.by %in% colnames(colData(sce)))
    col_data <- colData(sce)
    row_data <- rowData(sce)
    cells_split <- split(colnames(sce), colData(sce)[, split.by])
    # split assays
    sce_list <- lapply(cells_split, function(cells) {
        cells_index <- rownames(col_data) %in% cells
        inner_counts <- counts(sce)[, cells_index]
        inner_sce <- SingleCellExperiment(
            list(counts = inner_counts),
            rowData = row_data,
            colData = col_data[cells_index, ]
        )
        for (an in assayNames(sce)) {
            if (an == "counts") next()
            assay(inner_sce, an) <- assay(sce, an)[, cells_index]
        }
        inner_sce
    })
    sce_list
}
correct_batch <- function(sce, ...) {
    UseMethod("correct_batch", sce)
}
correct_batch.default <- function(sce, ...) {
    warning("No signatures for ", class(sce), ", ", "use seurat")
    class(sce) <- "seurat"
    correct_batch(sce, ...)
}
correct_batch.seurat <- function(sce, label = "batch", ...) {
    require(Seurat)
    class(sce) <- "SingleCellExperiment"
    infotag <- "[correct batch] "
    message(infotag, "create seurat object ...")
    scrna <- CreateSeuratObject(counts = counts(sce), project = "seurat", 
                                meta.data = data.frame(colData(sce)))
    scrna_split <- SplitObject(scrna, split.by = label)
    message(infotag, "do integration on each dataset ...")
    scrna_split <- lapply(scrna_split, function(sc) {
        NormalizeData(sc, verbose = FALSE) %>%
            ScaleData(verbose = FALSE) %>%
            FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    })
    scrna_anchors <- FindIntegrationAnchors(object.list = scrna_split, dims = 1:30)
    scrna_integrated <- IntegrateData(anchorset = scrna_anchors, dims = 1:30)
    DefaultAssay(scrna_integrated) <- "integrated"
    message(infotag, "do scale and pca on integrated dataset ...")
    scrna_integrated <- scrna_integrated %>% 
        ScaleData(verbose = F) %>%
        RunPCA(npcs = 30, verbose = FALSE) %>%
        RunTSNE(npcs = 1:20)
    scrna_integrated
}
correct_batch.harmony <- function(sce, label = "batch", ...) {
    require(harmony)
    class(sce) <- "SingleCellExperiment"
    infotag <- "[correct batch] "
    message(infotag, "get pcs from seurat3")
    require(Seurat)
    class(sce) <- "seurat"
    scrna_integrated <- correct_batch(sce, label, ...)
    ## below works as RunHarmony
    # scrna_pcs <- Embeddings(scrna_integrated, reduction = "pca")
    # harmony_pcs <- HarmonyMatrix(scrna_pcs, meta_data = scrna_integrated@meta.data,
    #                              vars_use = labname, do_pca = FALSE, npcs = 30,
    #                              max.iter.cluster = 50, max.iter.harnomy = 100)
    # Embeddings(scrna_integrated, "harmony") <- harmony_pcs
    scrna_integrated <- RunHarmony(scrna_integrated, label, assay.use = "integrated",
                                   max.iter.cluster = 50, max.iter.harnomy = 100)
    scrna_integrated <- RunTSNE(scrna_integrated, reduction = "harmony", dims = 1:20)
    scrna_integrated
}
correct_batch.liger <- function(sce, label = "batch", ...) {
    require(SingleCellExperiment)
    require(liger)
    class(sce) <- "SingleCellExperiment"
    infotag <- "[correct batch] "
    sce_split <- split_SCE(sce, label)
    counts_list <- lapply(sce_split, counts)
    message(infotag, "create, normalize and scale liger object ...")
    scrna <- createLiger(counts_list)
    scrna <- scrna %>% 
        normalize() %>% 
        selectGenes(do.plot = FALSE) %>%
        scaleNotCenter()
    message(infotag, "optimized and align datasets ...")
    # k_suggest <- suggestK(scrna, num.cores = ceiling(parallel::detectCores() / 8), 
    #                       gen.new = TRUE, return.results = TRUE ,plot.log2 = F, nrep = 5)
    scrna <- optimizeALS(scrna, k = 20, labmda = 5) %>%
        quantileAlignSNF(resolution = 0.4)
    scrna <- runTSNE(scrna, use.raw = FALSE)
    scrna
}
correct_batch.fastmnn <- function(sce, label = "batch", wrapper = TRUE, ...) {
    require(SingleCellExperiment)
    class(sce) <- "SingleCellExperiment"
    infotag <- "[correct batch] "
    message(infotag, "create seurat object ...")
    if (wrapper) {
        require(SeuratWrappers)
        require(Seurat)
        scrna <- CreateSeuratObject(counts = counts(sce), project = "seurat", 
                                    meta.data = data.frame(colData(sce)))
        scrna <- NormalizeData(scrna) %>% FindVariableFeatures()
        scrna_integrated <- SeuratWrappers::RunFastMNN(SplitObject(scrna, split.by = label))
        scrna_integrated <- RunTSNE(scrna_integrated, reduction = "mnn", dims = 1:30)
        scrna_integrated
    } else {
        require(scater)
        require(scran)
        message(infotag, " split and normalize batch dataset ...")
        sce_split <- split_SCE(sce, label)
        sce_split <- lapply(sce_split, scater::normalize)
        message(infotag, "find high variance genes")
        hvgs <- lapply(sce_split, function(sce) {
            # below two are used to find HVGs
            fit <- scran::trendVar(sce, parametric = TRUE, use.spikes = FALSE)
            decomp <- scran::decomposeVar(sce, fit)
            decomp
        })
        hvgs <- do.call(scran::combineVar, hvgs)
        top_hvgs <- order(hvgs$bio, decreasing = TRUE)[1:5000]
        logcounts_list <- lapply(sce_split, SingleCellExperiment::logcounts)
        message(infotag, "do cosine normalization....")
        cns <- lapply(logcounts_list, scran::cosineNorm)
        message(infotag, "do batch PCA")
        pcas <- do.call(scran::multiBatchPCA, c(cns, subset.row = top_hvgs, d = 50, BPPARAM = BiocSingular::RandomParam(deferred=TRUE)))
        message(infotag, "do fastmnn")
        mnnout <- do.call(scran::fastMNN, c(pcas, pc.input = TRUE, d = 50, BPPARAM = BiocSingular::RandomParam(deferred=TRUE)))# 30 componets
        mnnout
    }
}
correct_batch.scmerge <- function(sce, label = "batch", ...) {
    require(SingleCellExperiment)
    require(scater)
    require(scMerge)
    class(sce) <- "SingleCellExperiment"
    infotag <- "[correct batch] "
    sce_split <- split_SCE(sce, "batch")
    sce_split <- lapply(sce_split, scater::normalize)
    genes <- lapply(sce_split, function(inner_sce) {
        dat <- scSEGIndex(logcounts(inner_sce), ncore = 12)
        rownames(dat)
    })
    genes <- unlist(genes)
    sce <- scMerge::sce_cbind(sce_split)
    sce_merge <- scMerge(
        sce_combine = sce,
        kmeansK = c(11, 11),
        ctl = genes,
        assay_name = "scmerge_unsupervised",
        parallel = TRUE,
        fast_svd = TRUE)
    sce_merge <- scater::runPCA(sce_merge, ncomponents = 30, exprs_values = "scmerge_unsupervised")
    sce_merge
}
if (TRUE) {
    # correct_batch.zinbwave <- function(sce, label = "batch", ...) {
    #     require(zinbwave)
    #     require(SingleCellExperiment)
    #     require(Seurat)
    #     infotag <- "[correct batch] "
    #     message(infotag, "create seurat object ...")
    #     class(sce) <- "SingleCellExperiment"
    #     sce <- sce[rowSums(counts(sce) > 5) > 5, ]
    #     scrna <- CreateSeuratObject(counts = counts(sce), project = "seurat", 
    #                                 meta.data = data.frame(colData(sce)))
    #     scrna <- NormalizeData(scrna, verbose = FALSE) %>%
    #         ScaleData(verbose = FALSE) %>%
    #         FindVariableFeatures()
    #     sce <- SingleCellExperiment(
    #         list(counts = scrna@assays$RNA@counts,
    #              data = scrna@assays$RNA@data),
    #         colData = scrna@meta.data
    #     )
    #     sce <- zinbwave(sce, K = 2, epsilon=1000)
    #     # sce <- zinbsurf(sce, X = paste0("~ ", label), K = 2,
    #     #                 which_assay = "data", which_genes = VariableFeatures(scrna))
    #     sce
    # }
    # correct_batch.scanorama <- function(sce, label = "batch", ...) {
    #     require(SingleCellExperiment)
    #     require(reticulate)
    #     #reticulate::use_python("/home/xizhihui/biosoft/miniconda3/bin/python", required = TRUE)
    #     reticulate::use_virtualenv("/home/xizhihui/exercise/python_exercise/single_cell", required = TRUE)
    #     scanorama <- reticulate::import("scanorama")
    #     class(sce) <- "SingleCellExperiment"
    #     sce_split <- split_SCE(sce, label)
    #     datasets <- lapply(sce_split, counts)
    #     genes_list <- lapply(sce_split, rownames)
    #     # Integration and batch correction.
    #     scrna_integrated <- scanorama$correct(datasets, genes_list,
    #                                           return_dimred=TRUE, return_dense=TRUE)
    #     scrna_integrated
    # }
}

# load dataset2 ----
if (TRUE) {
    #https://raw.githubusercontent.com/JinmiaoChenLab/Batch-effect-removal-benchmarking/master/Data/dataset2/filtered_total_batch1_seqwell_batch2_10x.txt.gz
    dir.create("01_dataset", recursive = TRUE)
    dataset <- read.delim("filtered_total_batch1_seqwell_batch2_10x.txt")
    metadata <- stringr::str_split_fixed(colnames(dataset), "_", n = 3)
    metadata <- as.data.frame(metadata)
    metadata <- data.frame(celltype = metadata$V1,
                           batch = rep(c("seqwell", "10x"), c(4239, 2715))) #seq
    colnames(dataset) <- paste0(metadata$batch, 1:ncol(dataset))
    rownames(metadata) <- colnames(dataset)
    sce <- SingleCellExperiment(
        assays = list(counts = as.matrix(dataset)),
        colData = metadata
    )
    saveRDS(sce, file = "01_dataset/01_sce.rds")
}


# do batch correct ----
if (TRUE) {
    sce <- readRDS("01_dataset/01_sce.rds")
    sce_split <- split_SCE(sce, "batch")
    sce2 <- sce
    
    class(sce2) <- "seurat"
    sc_seurat <- correct_batch(sce2, label = "batch")
    saveRDS(sc_seurat, file = "01_dataset/02_sc_seurat.rds")
    class(sce2) <- "harmony"
    sc_harmony <- correct_batch(sce2, label = "batch")
    saveRDS(sc_harmony, file = "01_dataset/03_sc_harmony.rds")
    
    class(sce2) <- "liger"
    sc_liger <- correct_batch(sce2, label = "batch")
    saveRDS(sc_liger, file = "01_dataset/04_sc_liger.rds")
    
    class(sce2) <- "scmerge"
    sc_scmerge <- correct_batch(sce2, "batch")
    saveRDS(sc_scmerge, file = "01_dataset/05_sc_scmerge.rds")
    
    class(sce2) <- "fastmnn"
    sc_fastmnn <- correct_batch(sce2, "batch")
    saveRDS(sc_fastmnn, file = "01_dataset/06_sc_fastmnn.rds")
}

# do measure on correct results ----
if (TRUE) {
    get_pcas <- function(sce_list) {
        pcas <- list()
        methods <- names(sce_list)
        for (method in methods) {
            if (method %in% c("seurat", "harmony", "fastmnn")) {
                pca_name <- ifelse(method == "seurat", "pca",
                                   ifelse(method == "harmony", "harmony", "mnn"))
                pcas[[method]] <- Seurat::Embeddings(sce_list[[method]], pca_name)
            } else if (method == "liger") {
                pcas[[method]] <- sce_list[[method]]@H.norm
            } else if (method == "scmerge") {
                pcas[[method]] <- reducedDim(sce_list[[method]], "PCA")
            }
        }
        pcas
    }
    get_tsne <- function(sce_list, sce) {
        tsne_list <- list()
        methods <- names(sce_list)
        metadata <- as.data.frame(colData(sce))
        for (method in methods) {
            if (method %in% c("seurat", "harmony", "fastmnn")) {
                tsne <- Seurat::Embeddings(sce_list[[method]], "tsne")
            } else if (method == "liger") {
                tsne <- sce_list[[method]]@tsne.coords
            } else if (method == "scmerge") {
                pca <- reducedDim(sce_list[[method]], "PCA")
                tsne <- Rtsne::Rtsne(pca, num_threads = 12)
                tsne <- tsne$Y
                rownames(tsne) <- rownames(pca)
            }
            tsne <- as.data.frame(tsne)
            colnames(tsne) <- c("tSNE_1", "tSNE_2")
            tsne[, c("celltype", "batch")] <- metadata[match(rownames(tsne), rownames(metadata)), c("celltype", "batch")]
            tsne_list[[method]] <- tsne
        }
        tsne_list
    }
    do_kbet <- function(pcas) {
        kbet_list <- list()
        sce <- readRDS("01_dataset/01_sce.rds")
        methods <- names(pcas)
        for (method in methods) {
            message("start kBET on ", method)
            cur_pca <- pcas[[method]]
            cur_idx <- match(rownames(cur_pca), colnames(sce))
            cur_label <- colData(sce)[cur_idx, ]
            kbet_result <- apply(cur_label, 2, function(label) {
                result <- kBET::kBET(df = cur_pca[, 1:20], batch = label, do.pca = FALSE, plot = FALSE)
                result$average.pval
            })
            kbet_list[[method]] <- kbet_result
        }
        kbet_list
    }
    do_lisi <- function(pcas) {
        lisi_list <- list()
        sce <- readRDS("01_dataset/01_sce.rds")
        methods <- names(pcas)
        for (method in methods) {
            message("start LISI on ", method)
            cur_pca <- pcas[[method]]
            cur_idx <- match(rownames(cur_pca), colnames(sce))
            cur_label <- colData(sce)[cur_idx, ]
            result <- lisi::compute_lisi(cur_pca, cur_label, label_colnames = c("batch", "celltype"))
            lisi_list[[method]] <- apply(result, 2, function(x) {
                (median(x) - min(x)) / (max(x) - min(x))
            })
        }
        lisi_list
    }
    do_ari <- function(pcas) {
        ari_list <- list()
        sce <- readRDS("01_dataset/01_sce.rds")
        for (method in names(pcas)) {
            message("start ARI on ", method)
            cur_pca <- pcas[[method]]
            km <- kmeans(cur_pca, centers = 11)
            cluster <- km$cluster[match(rownames(cur_pca), names(km$cluster))]
            labels <- colData(sce)[match(rownames(cur_pca), colnames(sce)), ]
            batch_ari <- mclust::adjustedRandIndex(labels$batch, cluster)
            celltype_ari <- mclust::adjustedRandIndex(labels$celltype, cluster) 
            ari_list[[method]] = structure(c(batch_ari, celltype_ari), names = c("batch", "celltype"))
        }
        ari_list
    }
    do_asw <- function(pcas) {
        asw_list <- list()
        sce <- readRDS("01_dataset/01_sce.rds")
        methods <- names(pcas)
        for (method in methods) {
            message("start ASW on ", method)
            cur_pca <- pcas[[method]]
            cur_idx <- match(rownames(cur_pca), colnames(sce))
            cur_label <- colData(sce)[cur_idx, ]
            asw <- apply(cur_label, 2, function(cls) {
                cls <- as.numeric(as.factor(cls)) # need integer
                sk <- cluster::silhouette(cls, cluster::daisy(cur_pca))
                mean(sk[, 3])
            })
            asw_list[[method]] <- asw
        }
        asw_list
    }
    myscatter <- function(measures, mtype) {
        require(ggplot2)
        measurements <- as.data.frame(do.call(rbind, measures))
        measurements$method <- rownames(measurements)
        celltype_label <- paste(toupper(mtype), " cell type")
        batch_label <- paste0(toupper(mtype), " batch")
        if (mtype %in% c("asw", "ari", "kbet")) {
            measurements$batch <- 1 - measurements$batch
            batch_label <- paste0("1 - ", batch_label)
        } else if (mtype == "lisi") {
            measurements$celltype <- 1 - measurements$celltype
            celltype_label <- paste0("1 - ", celltype_label)
        }
        ggplot(measurements, aes(x = celltype, y = batch, color = method, label = method)) +
            geom_point() +
            ggrepel::geom_text_repel() +
            labs(x = celltype_label, y = batch_label) +
            theme(legend.position = "none")
    }
    mytsne <- function(tsne) {
        require(ggplot2)
        p1 <- ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, color = celltype)) + 
            geom_point(size = 0.1) +
            guides(colour = guide_legend(override.aes = list(size=4)))
        p2 <- ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, color = batch)) + 
            geom_point(size = 0.1) +
            guides(colour = guide_legend(override.aes = list(size=4)))
        cowplot::plot_grid(p1, p2)
    }
    
    sce <- readRDS("01_dataset/01_sce.rds")
    methods <- c("seurat", "harmony", "liger", "scmerge", "fastmnn")
    names(methods) <- methods
    sce_list <- lapply(methods, function(method) {
        filename <- paste0("0", which(method == methods) + 1, "_sc_", method, ".rds")
        print(filename)
        readRDS(file.path("01_dataset", filename))
    })
    
    pcas <- get_pcas(sce_list)
    measurements <- list(
        KBET = do_kbet(pcas),
        ARI = do_ari(pcas),
        ASW = do_asw(pcas),
        LISI = do_lisi(pcas)
    )
    myscatter(measurements$ARI , "ari")
    
    tsnes <- get_tsne(sce_list, sce)
    mytsne(tsnes$fastmnn)
}

