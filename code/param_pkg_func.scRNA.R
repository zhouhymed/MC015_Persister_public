library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggsci)
library(glue)
library(glmGamPoi)
library(openxlsx)
library(ggstatsplot)
library(viridis)

# Default parameters
# Seurat pre-processing
NormalizeData.METHOD <- "LogNormalize"
NormalizeData.ASSAY <- "RNA"
FindVariableFeatures.SELECT_METHOD <- "vst"
FindVariableFeatures.NFEATURES = 2000
RunPCA.NPCS <- 100
RunUMAP.DIMS <- 1:40
FindNeighbors.DIMS <- 1:40
FindClusters.RESOLUTION <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.3, 1.4, 1.5)
# FindMarkers
DE.MIN_CELL_NUM <- 15
DE.MIN_PCT.STRICT <- 0.25
DE.LFC_THRESHOLD.STRICT <- 1
DE.FDR_THRESHOLD.STRICT <- 0.05
DE.MIN_PCT.LOOSE <- 0.1
DE.LFC_THRESHOLD.LOOSE <- log2(1.5)
DE.FDR_THRESHOLD.LOOSE <- 0.1
# AddModuleScore
AddModuleScore.MINGSSIZE <- 3
# Plot
FeaturePlot.PAL_CONT <- viridis::viridis_pal(begin = 0.1, end = 0.9)(2)
# Over-representation enrichment
OR.PVALUE_THRESHOLD <- 0.05
OR.BAR.TOPN_GENESET <- 10
OR.SKIP_MSIGDB_GS_SUBCAT <- c("MIR:MIR_Legacy", "MIR:MIRDB", "GO:CC", "CM", "CGN", "HPO")
# tradeSeq
TRADESEQ.ASSOTEST.FDR_MAX <- 0.05
TRADESEQ.ASSOTEST.LFC_MIN <- log2(2)
TRADESEQ.SMOOTH.NPOINTS <- 100
TRADESEQ.SVE.FDR_MAX <- 0.05
TRADESEQ.SVE.LFC_MIN <- log2(1)
TRADESEQ.DIFFEND.LFC_MIN <- log2(1)
TRADESEQ.DIFFEND.FDR_MAX <- 0.05
TRADESEQ.PATTERN.LFC_MIN <- log2(1)
TRADESEQ.PATTERN.FDR_MAX <- 0.05
TRADESEQ.PATTERN.NPOINTS <- 100
TRADESEQ.PATTERN.FC_FILTER_MIN <- 0.5

#' Alias columns to metadata
add_alias_meta <- function(srt, mt_pattern = "^MT-", ribo_pattern = "^RP[SL]") {
  # srt$Project <- srt$orig.ident
  srt$nGene <- srt$nFeature_RNA
  srt$nUMI <- srt$nCount_RNA
  srt$log10GenesPerUMI <- log10(srt$nGene) / log10(srt$nUMI)
  srt$mitoRatio <-
    PercentageFeatureSet(object = srt, pattern = mt_pattern)
  srt$riboRatio <-
    PercentageFeatureSet(object = srt, pattern = ribo_pattern)
  srt
}


#' Standard preprocess for single sample
#'
#' Scale regress out cell cycle scores, mitoRatio and nFeature_RNA
standard_preprocess <-
  function(srt,
           NormalizeData.METHOD = "LogNormalize",
           NormalizeData.ASSAY = "RNA",
           FindVariableFeatures.SELECT_METHOD = "vst",
           FindVariableFeatures.NFEATURES = 2000,
           RunPCA.NPCS = 100,
           RunUMAP.DIMS = 1:40,
           FindNeighbors.DIMS = 1:40,
           FindClusters.RESOLUTION = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.3, 1.4, 1.5),
           prefix = NULL) {
    message("NormalizeData with default")
    srt <- NormalizeData(srt,
                         assay = NormalizeData.ASSAY,
                         normalization.method = NormalizeData.METHOD)
    gc()

    message("FindVariableFeatures with default")
    srt <-
      FindVariableFeatures(srt,
                           selection.method = FindVariableFeatures.SELECT_METHOD,
                           nfeatures = FindVariableFeatures.NFEATURES)
    gc()

    # Integrated assay need to keep values for cell cycle genes,
    # but standard preprocessing is okay.

    # Calcualte cell cycle scores
    message("CellCycleScoring with default")
    srt <-
      CellCycleScoring(srt,
                       g2m.features = cc.genes.updated.2019$g2m.genes,
                       s.features = cc.genes.updated.2019$s.genes)
    gc()
    # colnames(srt@meta.data)

    # Scale
    message("Scale: nFeature_RNA, mitoRatio, S.score, G2M.Score")
    srt <-
      ScaleData(srt,
                vars.to.regress = c("nFeature_RNA",
                                    "mitoRatio",
                                    "S.Score",
                                    "G2M.Score"))
    gc()

    # PCA
    message("RunPCA")
    DefaultAssay(srt) <- "RNA"
    srt <- RunPCA(srt, npcs = RunPCA.NPCS, verbose = FALSE)
    message(glue("Reduction: {paste(names(srt@reductions), collapse = ', ')}"))
    gc()

    # UMAP
    message("RunUMAP")
    srt <-
      RunUMAP(srt,
              dims = RunUMAP.DIMS,
              verbose = TRUE)
    message(glue("Reduction: {paste(names(srt@reductions), collapse = ', ')}"))
    gc()

    # PCA_SNN
    message("FindNeighbors, RNA.PCA_NN and RNA.PCA_SNN")
    DefaultAssay(srt) <- "RNA"
    srt <-
      FindNeighbors(srt,
                    reduction = "pca",
                    dims = FindNeighbors.DIMS,
                    graph.name = c("RNA.PCA_NN", "RNA.PCA_SNN"))
    message(glue("Graphs: {paste(names(srt@graphs), collapse = ', ')}"))
    gc()

    message("FindClusters from RNA.PCA_SNN")
    srt <-
      FindClusters(srt, resolution = FindClusters.RESOLUTION, graph.name = "RNA.PCA_SNN")
    message(glue("Clusters: {paste(standard_preprocess_resolution(srt), collapse = ', ')}"))
    used_res_name <- standard_preprocess_resolution(srt)
    gc()

    message("Update cluster label")
    for (reso in standard_preprocess_resolution(srt)) {
      srt@meta.data[, reso, drop = TRUE] %>%
        factor() %>%
        forcats::fct_relabel(., .fun = function(x) {
          paste0("C", as.integer(as.numeric(x) + 1))
        }) -> cluster_label
      # table(cluster_label)
      1:length(levels(cluster_label)) %>%
        paste0("C", .) -> new_levels
      cluster_label <- forcats::fct_relevel(cluster_label, new_levels)
      new_reso_name <- paste0(reso, ".label")
      srt[[new_reso_name]] <- cluster_label
    }
    message(glue("Clusters new label: {paste(grep(pattern = '.label$', colnames(srt@meta.data), value = TRUE), collapse = ', ')}"))

    srt_name <- paste0(prefix, ".srt")
    res_name <- paste0(prefix, ".res_name")
    assign(srt_name, srt)
    assign(res_name, used_res_name)

    srt
  }


#' Standard preprocess for single sample, resolution used
#'
#' Scale regress out cell cycle scores, mitoRatio and nFeature_RNA
standard_preprocess_resolution <- function(srt, pattern = "^RNA.PCA_SNN.res.[0-9.]+$") {
  grep(pattern = pattern, colnames(srt@meta.data), value = TRUE)
}

standard_preprocess_resolution_labelled <- function(srt, pattern = "^RNA.PCA_SNN.res.[0-9.]+.label$") {
  grep(pattern = pattern, colnames(srt@meta.data), value = TRUE)
}


#' Standard preprocess for several samples where cell cycle scores have been
#' calculated.
#'
#' Scale regress out cell cycle scores, mitoRatio and nFeature_RNA
standard_preprocess.merge <-
  function(srt,
           use_pre_cal_cell_cycle_scores = TRUE,
           prefix = NULL,
           Rda_file = NULL,
           NormalizeData.METHOD = "LogNormalize",
           NormalizeData.ASSAY = "RNA",
           FindVariableFeatures.SELECT_METHOD = "vst",
           FindVariableFeatures.NFEATURES = 2000,
           RunPCA.NPCS = 100,
           RunUMAP.DIMS = 1:40,
           FindNeighbors.DIMS = 1:40,
           FindClusters.RESOLUTION = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.3, 1.4, 1.5)) {
    message("NormalizeData with default")
    srt <- NormalizeData(srt,
                         assay = NormalizeData.ASSAY,
                         normalization.method = NormalizeData.METHOD)
    gc()

    message("FindVariableFeatures with default")
    srt <-
      FindVariableFeatures(srt,
                           selection.method = FindVariableFeatures.SELECT_METHOD,
                           nfeatures = FindVariableFeatures.NFEATURES)
    gc()

    # Integrated assay need to keep values for cell cycle genes,
    # but standard preprocessing is okay.

    # Calcualte cell cycle scores
    if (!use_pre_cal_cell_cycle_scores) {
      message("CellCycleScoring with default")
      srt <-
        CellCycleScoring(srt,
                         g2m.features = cc.genes.updated.2019$g2m.genes,
                         s.features = cc.genes.updated.2019$s.genes)
      gc()
    } else {
      message("Use pre-calculated cell cycle scores.")
    }

    # colnames(srt@meta.data)

    # Scale
    message("Scale: nFeature_RNA, mitoRatio, S.score, G2M.Score")
    srt <-
      ScaleData(srt,
                vars.to.regress = c("nFeature_RNA",
                                    "mitoRatio",
                                    "S.Score",
                                    "G2M.Score"))
    gc()

    # PCA
    message("RunPCA")
    DefaultAssay(srt) <- "RNA"
    srt <- RunPCA(srt, npcs = RunPCA.NPCS, verbose = FALSE)
    message(glue("Reduction: {paste(names(srt@reductions), collapse = ', ')}"))
    gc()

    # UMAP
    message("RunUMAP")
    srt <-
      RunUMAP(srt,
              dims = RunUMAP.DIMS,
              verbose = TRUE)
    message(glue("Reduction: {paste(names(srt@reductions), collapse = ', ')}"))
    gc()

    # PCA_SNN
    message("FindNeighbors, RNA.PCA_NN and RNA.PCA_SNN")
    DefaultAssay(srt) <- "RNA"
    srt <-
      FindNeighbors(srt,
                    reduction = "pca",
                    dims = FindNeighbors.DIMS,
                    graph.name = c("RNA.PCA_NN", "RNA.PCA_SNN"))
    message(glue("Graphs: {paste(names(srt@graphs), collapse = ', ')}"))
    gc()

    message("FindClusters from RNA.PCA_SNN")
    srt <-
      FindClusters(srt, resolution = FindClusters.RESOLUTION, graph.name = "RNA.PCA_SNN")
    message(glue("Clusters: {paste(standard_preprocess_resolution(srt), collapse = ', ')}"))
    used_res_name <- standard_preprocess_resolution(srt)
    gc()

    message("Update cluster label")
    for (reso in standard_preprocess_resolution(srt)) {
      srt@meta.data[, reso, drop = TRUE] %>%
        factor() %>%
        forcats::fct_relabel(., .fun = function(x) {
          paste0("C", as.integer(as.numeric(x) + 1))
        }) -> cluster_label
      # table(cluster_label)
      1:length(levels(cluster_label)) %>%
        paste0("C", .) -> new_levels
      cluster_label <- forcats::fct_relevel(cluster_label, new_levels)
      new_reso_name <- paste0(reso, ".label")
      srt[[new_reso_name]] <- cluster_label
    }
    message(glue("Clusters new label: {paste(grep(pattern = '.label$', colnames(srt@meta.data), value = TRUE), collapse = ', ')}"))

    srt_name <- paste0(prefix, ".srt")
    res_name <- paste0(prefix, ".res_name")
    assign(srt_name, srt)
    assign(res_name, used_res_name)

    # Save
    message("Save")
    save(
      list = c(
        srt_name,
        res_name,
        "NormalizeData.METHOD",
        "NormalizeData.ASSAY",
        "FindVariableFeatures.SELECT_METHOD",
        "FindVariableFeatures.NFEATURES",
        "RunPCA.NPCS",
        "RunUMAP.DIMS",
        "FindNeighbors.DIMS",
        "FindClusters.RESOLUTION"
      ),
      file = Rda_file
    )
    srt
  }


#' ggsci pal discrete
ggsci_pal_d <- function(ggsci_pal_name = "jco", pal_size = 5, ...) {
  # jco
  if (ggsci_pal_name == "jco") {
    if (pal_size > 10) {
      pal <- ggsci::pal_jco(...)(10)
      pal <- colorRampPalette(pal, ...)(pal_size)
    } else {
      pal <- ggsci::pal_jco(...)(pal_size)
    }
  } else if (ggsci_pal_name == "startrek") {
    # Star trek
    if (pal_size > 7) {
      pal <- ggsci::pal_startrek(...)(7)
      pal <- colorRampPalette(pal, ...)(pal_size)
    } else {
      pal <- ggsci::pal_startrek(...)(pal_size)
    }
  } else if (ggsci_pal_name == "aaas") {
    # aaas
    if (pal_size > 10) {
      pal <- ggsci::pal_aaas(...)(10)
      pal <- colorRampPalette(pal, ...)(pal_size)
    } else {
      pal <- ggsci::pal_aaas(...)(pal_size)
    }
  } else if (ggsci_pal_name == "npg") {
    # NPG
    if (pal_size > 10) {
      pal <- ggsci::pal_npg(...)(10)
      pal <- colorRampPalette(pal, ...)(pal_size)
    } else {
      pal <- ggsci::pal_npg(...)(pal_size)
    }
  } else {
    stop("pal not found")
  }
  pal
}
