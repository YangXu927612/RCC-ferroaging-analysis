################################################################################
# 19_scrnaseq_comprehensive_analysis.R
# Single-cell RNA-seq comprehensive analysis pipeline.
# - Quality control
# - Normalization and feature selection
# - PCA, clustering (UMAP)
# - Cell type annotation (manual + reference-based)
# - Marker gene analysis
# - Cell proportion analysis
# - Per-gene UMAP / heatmap / boxplot visualization
# - scRNAtoolVis downstream plots
# - In silico gene knockout via scTenifoldKnk
# All paths are RELATIVE; no absolute paths are used.
################################################################################
# Working directory (default = current working directory)
workDir <- getwd()
# ==============================================================================
# 01. Load Required Packages
# ==============================================================================
# Genes of interest (read from Intersect_Genes.csv in the working directory)
showGenes <- read.csv(file.path(workDir, "Intersect_Genes.csv"), stringsAsFactors = FALSE)[, 1]
showGenes <- trimws(showGenes)
showGenes <- showGenes[showGenes != ""]
# Define input file path (relative to the working directory)
input_file <- file.path(workDir, "scRNA_matrix.rds")

# Core single cell analysis packages / 核心单细胞分析包
library(Seurat)        # Single cell data analysis / 单细胞数据分析核心包
library(limma)         # Differential expression analysis / 差异表达分析
library(SingleR)       # Cell type annotation / 单细胞类型注释
library(celldex)       # Cell type reference datasets / 细胞类型注释数据集

library(DoubletFinder)
library(celda)          # decontX在celda包中
# Data processing packages / 数据处理包
library(dplyr)         # Data manipulation / 数据处理和管道操作
library(magrittr)      # Pipe operator / 管道操作符
library(tidyr)         # Data tidying / 数据整理
library(stringr)       # String processing / 字符串处理

# Visualization packages / 可视化包
library(ggplot2)       # Basic plotting / 基础绘图
library(ggpubr)        # Publication ready plots / 统计绘图和发表级图形
library(RColorBrewer)  # Color palettes / 调色板
library(viridis)       # Modern color schemes / 现代配色方案
library(scales)        # Scale functions / 图形标度
library(patchwork)     # Plot composition / 图形组合

# Statistical analysis packages / 统计分析包
library(rstatix)       # Statistical tests / 统计检验

# Interactive and utility packages / 交互和辅助包
library(DT)            # Interactive tables / 交互式表格
library(progress)      # Progress bars / 进度条显示

# File I/O packages / 文件输出包
library(openxlsx)      # Excel file I/O / Excel文件读写

# Gene clustering packages / 基因聚类包
library(ComplexHeatmap) # Complex heatmaps / 复杂热图
library(ggsci)         # Scientific journal color schemes / 科学期刊配色方案

# Data preprocessing related packages / 数据预处理相关包
library(Matrix)        # Sparse matrices / 稀疏矩阵
library(svglite)       # SVG export / 导出svg格式图
library(clustree)      # Clustering trees / 聚类树图
library(cluster)       # Silhouette coefficient / 轮廓系数计算
library(ggrepel)       # Avoid overlapping labels / 避免标签重叠
library(ggunchull)     # Convex hull for scRNAtoolVis / 凸包绘制
library(scRNAtoolVis)  # Single cell visualization tools / 单细胞可视化工具

# Advanced visualization packages / 高级可视化包
library(plot1cell)     # Circular UMAP plots / 环形UMAP图
library(circlize)      # Circular plot functions / 环形图绘制功能

message("All required R packages loaded successfully! / 所有必需的R包已成功加载完成！")

# ==============================================================================
# 02. Parameter Settings and Working Directory / 参数设置和工作目录
# ==============================================================================



# Check and set working directory / 检查并设置工作目录
if (!file.exists(workDir)) {
  stop("Working directory does not exist / 工作目录不存在，请检查路径：", workDir)
} else {
  setwd(workDir)
  message("Working directory set to / 工作目录设置为：", getwd())
}

# Analysis parameter settings / 分析参数设置
analysis_params <- list(
  # Basic parameters / 基础参数
  logFCfilter = 1,                    # logFC filtering threshold / logFC筛选阈值
  adjPvalFilter = 0.05,               # Adjusted p-value threshold / 调整后p值筛选阈值
  min_cells = 10,                     # Minimum cells expressing gene / 基因至少在多少个细胞中表达
  min_features = 30,                  # Minimum genes per cell / 细胞至少表达多少个基因

  # Quality control parameters / 质量控制参数
  nFeature_RNA_min = 20,              # Minimum genes per cell / 每个细胞最少基因数
  percent_mt_max = 20,                # Maximum mitochondrial percentage / 线粒体基因比例最大值(%)

  # Data processing parameters / 数据处理参数
  normalization_method = "LogNormalize", # Normalization method / 标准化方法
  scale_factor = 10000,               # Scale factor / 标准化缩放因子
  selection_method = "vst",           # Variable gene selection method / 高变基因选择方法
  n_variable_features = 1500,         # Number of variable genes / 高变基因数量

  # PCA and clustering parameters / PCA和聚类参数
  n_pcs = 20,                        # Number of principal components / PCA主成分数量
  pcSelect = 13,                     # PCs for clustering / 用于聚类的主成分数
  cluster_resolution = 0.2,          # Clustering resolution / 聚类分辨率
  cluster_resolution_range = seq(0.1, 1.0, 0.1), # Resolution range / 多分辨率聚类范围

  # Other parameters / 其他参数
  jackstraw_replicates = 100,        # JackStraw replicates / JackStraw重采样次数
  cumulative_variance_threshold = 0.90, # Cumulative variance threshold / 累积方差贡献率阈值
  small_cluster_threshold = 50,      # Small cluster threshold / 小聚类阈值
  max_clusters_silhouette = 20       # Max clusters for silhouette / 轮廓系数计算的最大聚类数
)



# Create output directory structure with numbering / 创建结果输出目录结构（带序号）
output_dirs <- list(
  data_preprocessing = "01.Data_Preprocessing",           # 数据预处理
  quality_control = "02.Quality_control",                # 质量控制
  feature_selection = "03.Feature_Selection",            # 特征选择
  pca_analysis = "04.PCA_Analysis",                      # PCA分析
  clustering = "05.Clustering_Analysis",                 # 聚类分析
  umap_visualization = "06.UMAP_Visualization",          # UMAP可视化
  cell_annotation = "07.Cell_Type_Annotation",           # 细胞类型注释
  marker_genes = "08.Marker_Gene_Analysis",              # 标记基因分析
  differential_analysis = "09.Differential_Gene_Analysis", # 差异基因分析
  proportion_analysis = "10.Cell_Proportion_Analysis",    # 细胞比例分析
  gene_expression = "11.Gene_Expression_Visualization",   # 基因表达可视化
  heatmaps = "12.Heatmap_Visualization",                 # 热图可视化
  statistics_reports = "13.Statistical_Reports",          # 统计报告
  boxplot_analysis = "14.BoxPlot_Analysis",              # 箱线图分析
  dotplot_analysis = "15.DotPlot_Analysis",              # 点图分析
  final_results = "16.Final_Results",                    # 最终结果
  individual_gene_umap = "17.Individual_Gene_UMAP",      # 单基因UMAP分析
  circular_umap = "19.Circular_UMAP",                     # 环形UMAP图
  scRNAtoolVis = "20.scRNAtoolVis_Visualization",          # scRNAtoolVis可视化
  gene_knockout = "21.Gene_Knockout_Analysis"              # 基因敲除模拟分析
)

# Create all output directories
message("Creating output directory structure...")
for (dir_name in output_dirs) {
  full_path <- file.path(workDir, dir_name)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
    message(sprintf("Created directory: %s", dir_name))
  }
}

# Create subdirectories
boxplot_subdirs <- c("ShowGenes_BoxPlots", "ShowGenes_Statistics", "ShowGenes_Individual_Plots")
for (subdir in boxplot_subdirs) {
  full_path <- file.path(workDir, output_dirs$boxplot_analysis, subdir)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
  }
}

# ==============================================================================
# 03. Data Reading and Processing
# ==============================================================================

message("Reading single cell data...")



# Check if file exists
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Read the RDS file directly
message(sprintf("Reading data from: %s", basename(input_file)))
single_cell_data <- readRDS(input_file)

# Check data dimensions
message(sprintf("Data loaded successfully: %d genes x %d cells",
                nrow(single_cell_data), ncol(single_cell_data)))

# Check and process group information
control_cells_found <- any(grepl("^control\\.", colnames(single_cell_data), ignore.case = FALSE))
disease_cells_found <- any(grepl("^Disease\\.", colnames(single_cell_data), ignore.case = FALSE))

if(control_cells_found) {
  message("Detected control prefix in cells")
}
if(disease_cells_found) {
  message("Detected Disease prefix in cells")
}

# Check if proper prefixes exist
if(!control_cells_found && !disease_cells_found) {
  stop("ERROR: No control or Disease prefixes found in cell names! Please ensure your data has cells with 'control.' or 'Disease.' prefixes")
}
 
# Verify both groups are present
if(!control_cells_found) {
  stop("ERROR: No control group cells found!")
}
if(!disease_cells_found) {
  stop("ERROR: No Disease group cells found!")
}

# Data quality check
if (ncol(single_cell_data) == 0) {
  stop("Data file is empty!")
}

message(sprintf("Data contains %d genes and %d cells", nrow(single_cell_data), ncol(single_cell_data)))
message(sprintf("数据包含 %d 个基因和 %d 个细胞", nrow(single_cell_data), ncol(single_cell_data)))

# Count group information
control_cells <- sum(grepl("^control\\.", colnames(single_cell_data)))
disease_cells <- sum(grepl("^Disease\\.", colnames(single_cell_data)))
message(sprintf("control cells: %d, Disease cells: %d", control_cells, disease_cells))
message(sprintf("对照组细胞数: %d, 疾病组细胞数: %d", control_cells, disease_cells))

# Save raw merged data
saveRDS(single_cell_data, file.path(output_dirs$data_preprocessing, "01_Raw_Combined_Data.rds"))

# ==============================================================================
# 04. Create Seurat Object and Quality control
# ==============================================================================

message("Creating Seurat object and performing quality control.")
message("创建Seurat对象并进行质量控制...")

# Create Seurat object
Peripheral_Blood_Mononuclear_Cells <- CreateSeuratObject(
  counts = single_cell_data,
  min.cells = analysis_params$min_cells,
  min.features = analysis_params$min_features
)

message("Seurat object created")

# Calculate mitochondrial percentage
Peripheral_Blood_Mononuclear_Cells[["percent.mt"]] <- PercentageFeatureSet(
  object = Peripheral_Blood_Mononuclear_Cells,
  pattern = "^MT-"
)

# Plot quality control violin plots BEFORE filtering
pdf(file = file.path(output_dirs$quality_control, "01_Quality_Violin_Before_Filter.pdf"),
    width = 18, height = 6.5)
print(
  VlnPlot(
    object = Peripheral_Blood_Mononuclear_Cells,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0
  )
)
dev.off()

# Statistics before filtering
cells_before <- ncol(Peripheral_Blood_Mononuclear_Cells)
genes_before <- nrow(Peripheral_Blood_Mononuclear_Cells)

# Data filtering
Peripheral_Blood_Mononuclear_Cells <- subset(
  x = Peripheral_Blood_Mononuclear_Cells,
  subset = nFeature_RNA > analysis_params$nFeature_RNA_min &
           percent.mt < analysis_params$percent_mt_max
)

# Plot quality control violin plots AFTER filtering
pdf(file = file.path(output_dirs$quality_control, "01_Quality_Violin_After_Filter.pdf"),
    width = 18, height = 6.5)
p_after <- VlnPlot(
  object = Peripheral_Blood_Mononuclear_Cells,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)
# Set y-axis limit for percent.mt plot to 100
p_after[[3]] <- p_after[[3]] + ylim(0, 100)
print(p_after)
dev.off()

# Statistics after filtering
cells_after <- ncol(Peripheral_Blood_Mononuclear_Cells)
genes_after <- nrow(Peripheral_Blood_Mononuclear_Cells)

message(sprintf("Before filtering: %d cells, %d genes", cells_before, genes_before))
message(sprintf("过滤前: %d 细胞, %d 基因", cells_before, genes_before))
message(sprintf("After filtering: %d cells, %d genes", cells_after, genes_after))
message(sprintf("过滤后: %d 细胞, %d 基因", cells_after, genes_after))
message(sprintf("Retention rate: cells %.2f%%, genes %.2f%%",
                (cells_after/cells_before)*100, (genes_after/genes_before)*100))
message(sprintf("保留率: 细胞 %.2f%%, 基因 %.2f%%",
                (cells_after/cells_before)*100, (genes_after/genes_before)*100))

# Save quality control summary
qc_summary <- data.frame(
  Metric = c("Cells_Before_Filter", "Cells_After_Filter", "Genes_Before_Filter",
             "Genes_After_Filter", "Cell_Retention_Rate", "Gene_Retention_Rate"),
  Value = c(cells_before, cells_after, genes_before, genes_after,
            paste0(round((cells_after/cells_before)*100, 2), "%"),
            paste0(round((genes_after/genes_before)*100, 2), "%"))
)

write.csv(qc_summary, file.path(output_dirs$quality_control, "02_QC_Summary.csv"),
          row.names = FALSE)

# ==============================================================================
# 04.3 Ambient RNA Removal (decontX)
# ==============================================================================

message("使用decontX去除环境RNA污染...")

# 获取原始count矩阵
raw_counts <- GetAssayData(Peripheral_Blood_Mononuclear_Cells, layer = "counts")

# 运行decontX
decontX_results <- decontX(raw_counts)

# 获取校正后的矩阵和污染比例
decontXcounts <- decontX_results$decontXcounts
contamination <- decontX_results$contamination

# 将污染比例添加到metadata
Peripheral_Blood_Mononuclear_Cells$decontX_contamination <- contamination

# 统计污染情况
message(sprintf("平均污染比例: %.2f%%", mean(contamination) * 100))
message(sprintf("中位污染比例: %.2f%%", median(contamination) * 100))
message(sprintf("最大污染比例: %.2f%%", max(contamination) * 100))

# 可视化污染分布
pdf(file.path(output_dirs$quality_control, "03_decontX_Contamination_Distribution.pdf"), width = 8, height = 6)
hist(contamination, breaks = 50, main = "decontX Contamination Distribution",
     xlab = "Contamination Proportion", col = "steelblue", border = "white")
abline(v = mean(contamination), col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Mean:", round(mean(contamination), 3)), col = "red", lty = 2, lwd = 2)
dev.off()

# 用校正后的矩阵替换原始数据
Peripheral_Blood_Mononuclear_Cells[["RNA"]]$counts <- decontXcounts

# 保存decontX结果摘要
decontX_summary <- data.frame(
  Metric = c("Mean_Contamination", "Median_Contamination", "Max_Contamination", "Min_Contamination"),
  Value = c(round(mean(contamination), 4), round(median(contamination), 4),
            round(max(contamination), 4), round(min(contamination), 4))
)
write.csv(decontX_summary, file.path(output_dirs$quality_control, "04_decontX_Summary.csv"), row.names = FALSE)

# 过滤高污染细胞 (污染比例 > 50%)
high_contamination_threshold <- 0.5
cells_before_contam_filter <- ncol(Peripheral_Blood_Mononuclear_Cells)

Peripheral_Blood_Mononuclear_Cells <- subset(
  Peripheral_Blood_Mononuclear_Cells,
  subset = decontX_contamination < high_contamination_threshold
)

cells_after_contam_filter <- ncol(Peripheral_Blood_Mononuclear_Cells)
message(sprintf("高污染细胞过滤: 去除 %d 个细胞 (污染>50%%)",
                cells_before_contam_filter - cells_after_contam_filter))
message(sprintf("过滤后剩余: %d 个细胞", cells_after_contam_filter))

message("decontX环境RNA去除完成！")

# ==============================================================================
# 04.5 Doublet Detection and Removal (DoubletFinder)
# 按分组分批运行，避免内存不足
# ==============================================================================

message("使用DoubletFinder进行双细胞检测和去除（分批处理）...")

# 记录去除前细胞数
cells_before_doublet <- ncol(Peripheral_Blood_Mononuclear_Cells)

# 添加分组信息用于分批
tmp_group <- ifelse(grepl("^control\\.", colnames(Peripheral_Blood_Mononuclear_Cells)), "control", "Disease")
Peripheral_Blood_Mononuclear_Cells$tmp_group <- tmp_group

# 分批处理函数
run_doubletfinder_batch <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

  # 参数扫描
  sweep.res <- paramSweep(seurat_obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_val <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  # 估算双细胞率
  n <- ncol(seurat_obj)
  dr <- n * 0.8 / 100000
  hp <- modelHomotypic(seurat_obj$seurat_clusters)
  nExp <- round(round(dr * n) * (1 - hp))

  # 运行DoubletFinder
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:10, pN = 0.25, pK = pK_val, nExp = nExp, sct = FALSE)

  df_col <- grep("^DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)[1]
  seurat_obj$Doublet_Classification <- seurat_obj@meta.data[[df_col]]
  return(seurat_obj)
}

# 分别对control和Disease运行
doublet_results <- list()
for (grp in c("control", "Disease")) {
  message(sprintf("  处理 %s 组...", grp))
  grp_obj <- subset(Peripheral_Blood_Mononuclear_Cells, subset = tmp_group == grp)
  message(sprintf("  %s 组细胞数: %d", grp, ncol(grp_obj)))

  tryCatch({
    grp_obj <- run_doubletfinder_batch(grp_obj)
    doublet_results[[grp]] <- grp_obj$Doublet_Classification
    names(doublet_results[[grp]]) <- colnames(grp_obj)
    n_doublet <- sum(grp_obj$Doublet_Classification == "Doublet")
    message(sprintf("  %s 组检测到 %d 个双细胞", grp, n_doublet))
  }, error = function(e) {
    stop(sprintf("%s 组DoubletFinder失败: %s", grp, e$message))
  })
}

# 合并结果
all_doublet_class <- c(doublet_results[["control"]], doublet_results[["Disease"]])
Peripheral_Blood_Mononuclear_Cells$Doublet_Classification <- all_doublet_class[colnames(Peripheral_Blood_Mononuclear_Cells)]

# 预处理用于UMAP可视化
Peripheral_Blood_Mononuclear_Cells <- NormalizeData(Peripheral_Blood_Mononuclear_Cells, verbose = FALSE)
Peripheral_Blood_Mononuclear_Cells <- FindVariableFeatures(Peripheral_Blood_Mononuclear_Cells, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
Peripheral_Blood_Mononuclear_Cells <- ScaleData(Peripheral_Blood_Mononuclear_Cells, verbose = FALSE)
Peripheral_Blood_Mononuclear_Cells <- RunPCA(Peripheral_Blood_Mononuclear_Cells, verbose = FALSE)
Peripheral_Blood_Mononuclear_Cells <- RunUMAP(Peripheral_Blood_Mononuclear_Cells, dims = 1:10, verbose = FALSE)

# 可视化双细胞分布
pdf(file.path(output_dirs$quality_control, "03_DoubletFinder_UMAP.pdf"), width = 10, height = 8)
print(DimPlot(Peripheral_Blood_Mononuclear_Cells, group.by = "Doublet_Classification", cols = c("Singlet" = "grey", "Doublet" = "red")))
dev.off()

# 去除双细胞
Peripheral_Blood_Mononuclear_Cells <- subset(Peripheral_Blood_Mononuclear_Cells, subset = Doublet_Classification == "Singlet")

# 清理临时列
Peripheral_Blood_Mononuclear_Cells$tmp_group <- NULL

# 统计结果
cells_after_doublet <- ncol(Peripheral_Blood_Mononuclear_Cells)
message(sprintf("双细胞去除前: %d 细胞", cells_before_doublet))
message(sprintf("双细胞去除后: %d 细胞", cells_after_doublet))
message(sprintf("去除双细胞: %d (%.2f%%)", cells_before_doublet - cells_after_doublet, (cells_before_doublet - cells_after_doublet)/cells_before_doublet*100))

# 重置对象用于后续标准分析
Peripheral_Blood_Mononuclear_Cells@reductions <- list()
Peripheral_Blood_Mononuclear_Cells@graphs <- list()
Peripheral_Blood_Mononuclear_Cells$seurat_clusters <- NULL

# ==============================================================================
# 05. Data Preprocessing
# ==============================================================================

message("Starting data preprocessing")

# Data normalization
Peripheral_Blood_Mononuclear_Cells <- NormalizeData(
  object = Peripheral_Blood_Mononuclear_Cells,
  normalization.method = analysis_params$normalization_method,
  scale.factor = analysis_params$scale_factor
)
message("Data normalization completed")

# Find variable features
Peripheral_Blood_Mononuclear_Cells <- FindVariableFeatures(
  object = Peripheral_Blood_Mononuclear_Cells,
  selection.method = analysis_params$selection_method,
  nfeatures = analysis_params$n_variable_features
)
message(sprintf("Identified %d variable features", length(VariableFeatures(Peripheral_Blood_Mononuclear_Cells))))
message(sprintf("识别到 %d 个高变基因", length(VariableFeatures(Peripheral_Blood_Mononuclear_Cells))))

# View top 10 variable genes
top10_genes <- head(VariableFeatures(Peripheral_Blood_Mononuclear_Cells), 10)
message("Top 10 variable genes")

# Plot feature selection
pdf(file = file.path(output_dirs$feature_selection, "01_Variable_Features.pdf"),
    width = 15, height = 8)
plot1 <- VariableFeaturePlot(Peripheral_Blood_Mononuclear_Cells)
plot2 <- LabelPoints(plot = plot1, points = top10_genes, repel = TRUE)
print(plot1 + plot2)
dev.off()

# Save variable genes list
write.csv(
  data.frame(Gene = VariableFeatures(Peripheral_Blood_Mononuclear_Cells)),
  file.path(output_dirs$feature_selection, "02_Variable_Genes_List.csv"),
  row.names = FALSE
)

# Scale data
Peripheral_Blood_Mononuclear_Cells <- ScaleData(Peripheral_Blood_Mononuclear_Cells)
message("Data scaling completed")

# ==============================================================================
# 06. PCA Analysis
# ==============================================================================

message("Performing PCA analysis")

# Run PCA
Peripheral_Blood_Mononuclear_Cells <- RunPCA(
  object = Peripheral_Blood_Mononuclear_Cells,
  npcs = analysis_params$n_pcs,
  features = VariableFeatures(object = Peripheral_Blood_Mononuclear_Cells)
)

# PCA visualization
pdf(file = file.path(output_dirs$pca_analysis, "01_PCA_Plot.pdf"),
    width = 8, height = 6)
print(DimPlot(Peripheral_Blood_Mononuclear_Cells, reduction = "pca") +
      NoGrid())  # Remove grid lines but keep axes
dev.off()

# PCA heatmaps (first 4 PCs)
for (i in 1:4) {
  pdf(file = file.path(output_dirs$pca_analysis, paste0("02_PCA_Heatmap_PC", i, ".pdf")),
      width = 10, height = 8)
  print(
    DimHeatmap(
      Peripheral_Blood_Mononuclear_Cells,
      dims = i,
      cells = 500,
      balanced = TRUE
    )
  )
  dev.off()
}

# JackStraw analysis
Peripheral_Blood_Mononuclear_Cells <- JackStraw(
  Peripheral_Blood_Mononuclear_Cells,
  num.replicate = analysis_params$jackstraw_replicates
)
Peripheral_Blood_Mononuclear_Cells <- ScoreJackStraw(
  Peripheral_Blood_Mononuclear_Cells,
  dims = 1:analysis_params$n_pcs
)

# JackStraw plot
pdf(file = file.path(output_dirs$pca_analysis, "03_JackStraw_Plot.pdf"),
    width = 9, height = 6)
print(JackStrawPlot(Peripheral_Blood_Mononuclear_Cells, dims = 1:15))
dev.off()

# Save PCA results
pca_results <- Peripheral_Blood_Mononuclear_Cells[["pca"]]@stdev
pca_variance <- pca_results^2
cumulative_variance <- cumsum(pca_variance)

pca_table <- data.frame(
  PC = 1:length(pca_variance),
  StdDev = pca_results,
  Variance = pca_variance,
  CumulativeVariance = cumulative_variance
)

write.csv(pca_table, file.path(output_dirs$pca_analysis, "04_PCA_Results.csv"),
          row.names = FALSE)

# Select principal components
selected_pcs <- which(cumulative_variance >= analysis_params$cumulative_variance_threshold)[1]
if(is.na(selected_pcs)) selected_pcs <- analysis_params$pcSelect
message(sprintf("Selected %d principal components based on cumulative variance", selected_pcs))
message(sprintf("根据累积方差选择 %d 个主成分", selected_pcs))

# ==============================================================================
# 07. Clustering Analysis
# ==============================================================================

message("Performing clustering analysis")

# Find neighbors
Peripheral_Blood_Mononuclear_Cells <- FindNeighbors(
  object = Peripheral_Blood_Mononuclear_Cells,
  dims = 1:analysis_params$pcSelect
)

# Multi-resolution clustering
message("Performing multi-resolution clustering")
for (res in analysis_params$cluster_resolution_range) {
  Peripheral_Blood_Mononuclear_Cells <- FindClusters(
    Peripheral_Blood_Mononuclear_Cells,
    resolution = res,
    verbose = FALSE
  )
}

# Generate Clustree plot
pdf(file.path(output_dirs$clustering, "01_Clustree_Analysis.pdf"),
    width = 15, height = 10)
clustree_plot <- clustree(Peripheral_Blood_Mononuclear_Cells@meta.data, prefix = "RNA_snn_res.") +
  theme_classic() +
  labs(title = "Clustering Tree Across Different Resolutions")
print(clustree_plot)
dev.off()

# Use default resolution for final clustering
Peripheral_Blood_Mononuclear_Cells <- FindClusters(
  object = Peripheral_Blood_Mononuclear_Cells,
  resolution = analysis_params$cluster_resolution
)
message(sprintf("Clustering completed with resolution %.1f, identified %d clusters",
                analysis_params$cluster_resolution,
                length(unique(Peripheral_Blood_Mononuclear_Cells$seurat_clusters))))
message(sprintf("使用分辨率 %.1f 完成聚类，识别到 %d 个聚类",
                analysis_params$cluster_resolution,
                length(unique(Peripheral_Blood_Mononuclear_Cells$seurat_clusters))))

# UMAP dimension reduction
Peripheral_Blood_Mononuclear_Cells <- RunUMAP(
  object = Peripheral_Blood_Mononuclear_Cells,
  dims = 1:analysis_params$pcSelect
)

# ==============================================================================
# 08. UMAP Visualization
# ==============================================================================

message("Generating UMAP visualization")

# UMAP clustering plot
pdf(file = file.path(output_dirs$umap_visualization, "01_UMAP_Clusters.pdf"),
    width = 10, height = 8)
print(
  DimPlot(
    object = Peripheral_Blood_Mononuclear_Cells,
    reduction = "umap",
    pt.size = 1.5,
    label = TRUE,
    label.size = 4
  ) +
  ggtitle("UMAP Clustering Results") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  )
)
dev.off()

# Save clustering information
write.table(
  Peripheral_Blood_Mononuclear_Cells$seurat_clusters,
  file = file.path(output_dirs$clustering, "02_Cluster_Assignments.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE
)

# ==============================================================================
# 09. Automatic Cell Type Annotation (SingleR)
# ==============================================================================

message("Performing automatic cell type annotation using SingleR...")
message("使用SingleR进行自动细胞类型注释...")

# Check if SingleR and celldex are available
if (!requireNamespace("SingleR", quietly = TRUE) || !requireNamespace("celldex", quietly = TRUE)) {
  message("Warning: SingleR or celldex package not installed, using default annotation")
  message("警告：SingleR 或 celldex 包未安装，使用默认注释")

  # Create default annotations
  default_labels <- paste0("Cluster_", Peripheral_Blood_Mononuclear_Cells$seurat_clusters)
  Peripheral_Blood_Mononuclear_Cells$SingleR_labels <- default_labels
  Peripheral_Blood_Mononuclear_Cells$cell_type <- default_labels
  Idents(Peripheral_Blood_Mononuclear_Cells) <- default_labels

} else {
  # CORRECT SingleR annotation workflow based on reference code
  # 基于参考代码的正确SingleR注释流程

  tryCatch({
    # Step 1: Extract expression data and cluster information
    # 步骤1：提取表达数据和聚类信息
    message("Extracting expression data and cluster information...")
    message("提取表达数据和聚类信息...")

    counts <- GetAssayData(object = Peripheral_Blood_Mononuclear_Cells, layer = "data")
    clusters <- Peripheral_Blood_Mononuclear_Cells@meta.data$seurat_clusters

    message("Expression matrix dimensions: ", nrow(counts), " genes x ", ncol(counts), " cells")
    message("表达矩阵维度: ", nrow(counts), " 基因 x ", ncol(counts), " 细胞")
    message("Number of clusters: ", length(unique(clusters)))
    message("聚类数量: ", length(unique(clusters)))

    # Step 2: Load reference dataset
    # 步骤2：加载参考数据集
    message("Loading HumanPrimaryCellAtlasData reference...")
    message("加载HumanPrimaryCellAtlasData参考数据集...")
    ref <- celldex::HumanPrimaryCellAtlasData()

    message("Reference dataset loaded: ", nrow(ref), " genes x ", ncol(ref), " samples")
    message("参考数据集已加载: ", nrow(ref), " 基因 x ", ncol(ref), " 样本")

    # Step 3: Perform cluster-level annotation (more robust)
    # 步骤3：进行聚类水平注释（更稳健）
    message("Running SingleR cluster-level annotation...")
    message("运行SingleR聚类水平注释...")

    singler_clusters <- SingleR(
      test = counts,
      ref = ref,
      labels = ref$label.main,
      clusters = clusters
    )

    # Save cluster annotation results
    # 保存聚类注释结果
    clusterAnn <- data.frame(
      Cluster = rownames(singler_clusters),
      CellType = singler_clusters$labels,
      stringsAsFactors = FALSE
    )

    write.csv(clusterAnn,
              file = file.path(output_dirs$cell_annotation, "Cluster_to_CellType_Mapping.csv"),
              row.names = FALSE)

    message("Cluster annotation completed. Results:")
    message("聚类注释完成。结果:")
    for(i in 1:nrow(clusterAnn)) {
      message(sprintf("  Cluster %s -> %s", clusterAnn$Cluster[i], clusterAnn$CellType[i]))
    }

    # Step 4: Map cluster annotations to individual cells
    # 步骤4：将聚类注释映射到单个细胞
    message("Mapping cluster annotations to individual cells...")
    message("将聚类注释映射到单个细胞...")

    # Create cluster to cell type mapping
    # 创建聚类到细胞类型的映射
    cluster_to_celltype <- setNames(clusterAnn$CellType, clusterAnn$Cluster)

    # Get cluster information for each cell
    # 获取每个细胞的聚类信息
    cell_clusters <- Peripheral_Blood_Mononuclear_Cells@meta.data$seurat_clusters

    # Map cluster annotations to each cell
    # 将聚类注释映射到每个细胞
    cell_types_from_clusters <- cluster_to_celltype[as.character(cell_clusters)]
    names(cell_types_from_clusters) <- colnames(Peripheral_Blood_Mononuclear_Cells)

    # Check for any missing annotations
    if(any(is.na(cell_types_from_clusters))) {
      message("Warning: Some cells have missing annotations, using cluster names as fallback")
      message("警告：部分细胞缺失注释，使用聚类名称作为备选")
      cell_types_from_clusters[is.na(cell_types_from_clusters)] <- paste0("Cluster_",
        cell_clusters[is.na(cell_types_from_clusters)])
    }

    # Step 5: Add annotations to Seurat object
    # 步骤5：将注释添加到Seurat对象
    message("Adding cell type annotations to Seurat object...")
    message("将细胞类型注释添加到Seurat对象...")

    # Add cluster-based cell type annotation to metadata
    # 将基于聚类的细胞类型注释添加到元数据中
    Peripheral_Blood_Mononuclear_Cells <- AddMetaData(
      Peripheral_Blood_Mononuclear_Cells,
      metadata = cell_types_from_clusters,
      col.name = "cell_type"
    )

    # Also add as SingleR_labels for compatibility
    # 同时添加为SingleR_labels以保持兼容性
    Peripheral_Blood_Mononuclear_Cells$SingleR_labels <- cell_types_from_clusters

    # Set Idents to cell types for downstream analysis
    # 设置Idents为细胞类型用于下游分析
    Idents(Peripheral_Blood_Mononuclear_Cells) <- cell_types_from_clusters

    # Verify annotations were added successfully
    # 验证注释是否成功添加
    message("Verification:")
    message("验证:")
    message("  SingleR_labels column exists: ", "SingleR_labels" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data))
    message("  cell_type column exists: ", "cell_type" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data))
    message("  Number of unique cell types: ", length(unique(cell_types_from_clusters)))
    message("  唯一细胞类型数量: ", length(unique(cell_types_from_clusters)))

    # Print cell type statistics
    # 打印细胞类型统计
    cell_type_table <- table(cell_types_from_clusters)
    message("\nCell type distribution:")
    message("细胞类型分布:")
    for(ct_name in names(cell_type_table)) {
      message(sprintf("  %s: %d cells", ct_name, cell_type_table[ct_name]))
    }

    message("SingleR annotation completed successfully!")
    message("SingleR注释成功完成！")

  }, error = function(e) {
    message("ERROR: SingleR annotation failed, using default annotation")
    message(sprintf("错误：SingleR注释失败，使用默认注释。错误信息: %s", e$message))
    message("Error details: ", conditionMessage(e))
    message("错误详情: ", conditionMessage(e))

    # Create default annotations on error
    message("Creating default cluster-based annotations...")
    message("创建默认的基于聚类的注释...")

    default_labels <- paste0("Cluster_", Peripheral_Blood_Mononuclear_Cells$seurat_clusters)
    Peripheral_Blood_Mononuclear_Cells$SingleR_labels <- default_labels
    Peripheral_Blood_Mononuclear_Cells$cell_type <- default_labels
    Idents(Peripheral_Blood_Mononuclear_Cells) <- default_labels

    # Verify default annotations were added
    message("Default annotation verification:")
    message("默认注释验证:")
    message("  SingleR_labels column exists: ", "SingleR_labels" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data))
    message("  cell_type column exists: ", "cell_type" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data))
  })
}

# ==============================================================================
# 09.5 Apply Manual Cell Type Annotations (优先于 SingleR 自动注释)
# 09.5 应用手动细胞类型注释（覆盖 SingleR 自动注释）
# ==============================================================================

message("\n================================================================================")
message("         Step 09.5: Applying Manual Cell Type Annotations                       ")
message("         步骤 09.5: 应用手动细胞类型注释（覆盖自动注释）                          ")
message("================================================================================")

# ------------------------------------------------------------------
# 修复: 提前创建 Type 列 (基于细胞名 control./Disease. 前缀)
# 原始代码在第963行才添加 Type, 但手动注释代码需要用到它
# ------------------------------------------------------------------
if (!"Type" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data)) {
  Type_tmp <- ifelse(
    grepl("^control\\.", colnames(Peripheral_Blood_Mononuclear_Cells)),
    "control", "Disease"
  )
  names(Type_tmp) <- colnames(Peripheral_Blood_Mononuclear_Cells)
  Peripheral_Blood_Mononuclear_Cells <- AddMetaData(
    object = Peripheral_Blood_Mononuclear_Cells,
    metadata = Type_tmp,
    col.name = "Type"
  )
  message("  ✓ 已提前创建 Type 列 (control/Disease)")
  message("  ✓ Type column created in advance")
}

# 查找手动注释文件 (优先按以下顺序: 细胞级 > 聚类级)
manual_annotation_files <- c(
  "Manual_Cell_Type_Annotations_7_CellTypes.csv",  # 细胞级注释
  "Manual_Cluster_to_7_CellTypes_Mapping.csv"      # 聚类级映射
)

manual_annotation_file <- NULL
manual_annotation_type <- NULL
for (f in manual_annotation_files) {
  fp <- file.path(output_dirs$cell_annotation, f)
  if (file.exists(fp)) {
    manual_annotation_file <- fp
    manual_annotation_type <- ifelse(grepl("Cluster_to", f), "cluster", "cell")
    message(sprintf("  找到手动注释文件: %s [%s级]", f, manual_annotation_type))
    message(sprintf("  Found manual annotation file: %s [%s-level]", f, manual_annotation_type))
    break
  }
}

if (!is.null(manual_annotation_file)) {
  manual_data <- read.csv(manual_annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
  message(sprintf("  加载手动注释: %d 行, %d 列", nrow(manual_data), ncol(manual_data)))
  message(sprintf("  列名: %s", paste(colnames(manual_data), collapse = ", ")))

  manual_cell_types <- NULL

  # Case 1: 细胞级注释文件 (有 Cell_ID + Manual_celltype_7)
  if ("Cell_ID" %in% colnames(manual_data) && "Manual_celltype_7" %in% colnames(manual_data)) {
    message("  -> 使用细胞级手动注释 (Cell-level manual annotation)")

    # 去除可能由 write.csv 添加的引号
    manual_data$Cell_ID <- as.character(manual_data$Cell_ID)
    manual_data$Cell_ID <- gsub('^"|"$', '', manual_data$Cell_ID)

    # 通过 Cell_ID 匹配
    cells_in_data <- colnames(Peripheral_Blood_Mononuclear_Cells)
    cells_match <- intersect(manual_data$Cell_ID, cells_in_data)

    message(sprintf("  手动注释包含 %d 个细胞ID, 当前数据有 %d 个细胞",
                    nrow(manual_data), length(cells_in_data)))
    message(sprintf("  Cell ID匹配: %d / %d", length(cells_match), length(cells_in_data)))

    if (length(cells_match) > 0) {
      cell_type_lookup <- setNames(manual_data$Manual_celltype_7, manual_data$Cell_ID)
      manual_cell_types <- cell_type_lookup[cells_in_data]
      names(manual_cell_types) <- cells_in_data

      # 去除 NA
      manual_cell_types <- as.character(manual_cell_types)
      n_na_manual <- sum(is.na(manual_cell_types))
      if (n_na_manual > 0) {
        message(sprintf("  警告: %d 个细胞在手动注释中缺失, 将使用 SingleR 注释填充",
                        n_na_manual))
      }
    } else {
      message("  ⚠ Cell_ID 无法直接匹配, 回退到基于聚类匹配")
    }
  }

  # Case 2: 聚类级映射文件 (有 Cluster + Manual_celltype_7)
  if (is.null(manual_cell_types) &&
      "Cluster" %in% colnames(manual_data) &&
      "Manual_celltype_7" %in% colnames(manual_data)) {
    message("  -> 使用聚类级手动注释 (Cluster-level manual annotation)")

    # 标准化聚类ID (去除可能的引号)
    manual_data$Cluster <- as.character(manual_data$Cluster)
    manual_data$Cluster <- gsub('^"|"$', '', manual_data$Cluster)

    # 创建聚类到细胞类型的命名向量
    cluster_to_celltype <- setNames(
      manual_data$Manual_celltype_7,
      manual_data$Cluster
    )

    # 应用到每个细胞
    cell_clusters <- as.character(Peripheral_Blood_Mononuclear_Cells$seurat_clusters)
    manual_cell_types <- cluster_to_celltype[cell_clusters]
    names(manual_cell_types) <- colnames(Peripheral_Blood_Mononuclear_Cells)

    message(sprintf("  应用聚类级注释到 %d 个细胞", length(manual_cell_types)))
    message(sprintf("  聚类映射规则数: %d", length(cluster_to_celltype)))

    # 调试: 输出聚类映射
    for (cl in sort(unique(manual_data$Cluster))) {
      ct <- unique(manual_data$Manual_celltype_7[manual_data$Cluster == cl])
      message(sprintf("    Cluster %s -> %s", cl, ct))
    }
  }

  # 应用手动注释到 Seurat 对象
  if (!is.null(manual_cell_types)) {
    manual_cell_types <- as.character(manual_cell_types)

    # 处理 NA: 用 SingleR 结果填充
    n_na <- sum(is.na(manual_cell_types))
    if (n_na > 0) {
      singler_fallback <- as.character(Peripheral_Blood_Mononuclear_Cells$SingleR_labels)
      manual_cell_types[is.na(manual_cell_types)] <- singler_fallback[is.na(manual_cell_types)]
      message(sprintf("  已用 SingleR 注释填充 %d 个缺失的细胞类型", n_na))
    }

    # 添加到 metadata
    Peripheral_Blood_Mononuclear_Cells$Manual_celltype_7 <- manual_cell_types
    Peripheral_Blood_Mononuclear_Cells$cell_type <- manual_cell_types
    Peripheral_Blood_Mononuclear_Cells$SingleR_labels <- manual_cell_types

    # 设置 Idents 为手动注释的细胞类型
    Idents(Peripheral_Blood_Mononuclear_Cells) <- manual_cell_types

    # 输出统计
    manual_counts <- table(manual_cell_types)
    message("\n  手动注释细胞类型分布:")
    message("  Manual annotation cell type distribution:")
    for (ct in names(sort(manual_counts, decreasing = TRUE))) {
      message(sprintf("    %s: %d cells (%.2f%%)",
                      ct, manual_counts[ct],
                      manual_counts[ct] / sum(manual_counts) * 100))
    }

    # 保存完整手动注释结果
    manual_ann_out <- data.frame(
      Cell_ID = colnames(Peripheral_Blood_Mononuclear_Cells),
      Cluster = as.character(Peripheral_Blood_Mononuclear_Cells$seurat_clusters),
      SingleR_Label = as.character(Peripheral_Blood_Mononuclear_Cells$SingleR_labels),
      Manual_celltype_7 = manual_cell_types,
      Group = Peripheral_Blood_Mononuclear_Cells$Type,
      stringsAsFactors = FALSE
    )

    write.csv(
      manual_ann_out,
      file.path(output_dirs$cell_annotation, "05_Manual_Annotation_Applied.csv"),
      row.names = FALSE
    )

    # 重新绘制 UMAP (基于手动注释)
    tryCatch({
      pdf(file = file.path(output_dirs$cell_annotation, "06_Manual_Annotation_UMAP.pdf"),
          width = 12, height = 8)
      print(
        DimPlot(
          Peripheral_Blood_Mononuclear_Cells,
          reduction = "umap",
          pt.size = 1.5,
          label = TRUE,
          repel = TRUE,
          group.by = "Manual_celltype_7"
        ) +
          ggtitle("Manual Cell Type Annotation (7 Cell Types)") +
          theme_minimal() +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"),
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black", face = "bold")
          )
      )
      dev.off()

      # 分组UMAP
      pdf(file = file.path(output_dirs$cell_annotation, "07_Manual_Annotation_UMAP_Split.pdf"),
          width = 16, height = 7)
      print(
        DimPlot(
          Peripheral_Blood_Mononuclear_Cells,
          reduction = "umap",
          pt.size = 1,
          label = TRUE,
          repel = TRUE,
          group.by = "Manual_celltype_7",
          split.by = "Type"
        ) +
          ggtitle("Manual Cell Type Annotation: control vs Disease") +
          theme_minimal() +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"),
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black", face = "bold")
          )
      )
      dev.off()
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message(sprintf("手动注释UMAP绘图失败: %s", e$message))
    })

    message("  ✓ 手动注释已成功应用并保存!")
    message("  ✓ Manual annotation applied successfully!")
  } else {
    message("  ✗ 手动注释无法解析，保持 SingleR 自动注释")
    message("  ✗ Failed to parse manual annotation, retaining SingleR")
  }
} else {
  message("  未找到手动注释文件 (Manual_*.csv), 使用 SingleR 自动注释")
  message("  No manual annotation file found, using SingleR auto-annotation")
}

message("================================================================================\n")

message(sprintf("Automatic annotation completed, identified %d cell types",
                length(unique(Idents(Peripheral_Blood_Mononuclear_Cells)))))
message(sprintf("自动注释完成，识别到 %d 种细胞类型",
                length(unique(Idents(Peripheral_Blood_Mononuclear_Cells)))))

# Print cell type statistics
cell_type_counts <- table(Idents(Peripheral_Blood_Mononuclear_Cells))
message("\nCell type statistics")
for(ct in names(sort(cell_type_counts, decreasing = TRUE))) {
  message(sprintf("  %s: %d cells", ct, cell_type_counts[ct]))
}

# Visualize annotation results
pdf(file = file.path(output_dirs$cell_annotation, "01_Auto_Annotation_UMAP.pdf"),
    width = 12, height = 8)
print(
  DimPlot(
    Peripheral_Blood_Mononuclear_Cells,
    reduction = "umap",
    pt.size = 1.5,
    label = TRUE,
    repel = TRUE
  ) +
  ggtitle("Automatic Cell Type Annotation (SingleR)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  )
)
dev.off()

# Save annotation results
# ROBUST ANNOTATION HANDLING - Multiple layers of fallback
message("Starting robust annotation result processing...")
message("开始强健的注释结果处理...")

# Step 1: Ensure cell_type column exists
if(!"cell_type" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data)) {
  message("cell_type column not found, creating from current Idents...")
  message("未找到cell_type列，从当前Idents创建...")

  # Force create the cell_type column from Idents
  tryCatch({
    Peripheral_Blood_Mononuclear_Cells$cell_type <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
    message("Successfully created cell_type column from Idents")
    message("成功从Idents创建cell_type列")
  }, error = function(e) {
    message("Error creating cell_type from Idents: ", e$message)
    message("从Idents创建cell_type时出错: ", e$message)
    # Create from cluster numbers as absolute fallback
    Peripheral_Blood_Mononuclear_Cells$cell_type <<- paste0("Cluster_", Peripheral_Blood_Mononuclear_Cells$seurat_clusters)
  })
}

# Step 2: Verify and extract cell type data with multiple fallbacks
cell_type_data <- NULL

# First priority: cell_type column
if("cell_type" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data)) {
  tryCatch({
    cell_type_data <- as.character(Peripheral_Blood_Mononuclear_Cells$cell_type)
    message("Successfully extracted cell_type data from cell_type column")
    message("成功从cell_type列提取细胞类型数据")
  }, error = function(e) {
    message("Error accessing cell_type column: ", e$message)
    cell_type_data <- NULL
  })
}

# Second priority: SingleR_labels column
if(is.null(cell_type_data) && "SingleR_labels" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data)) {
  tryCatch({
    cell_type_data <- as.character(Peripheral_Blood_Mononuclear_Cells$SingleR_labels)
    message("Using SingleR_labels as cell_type")
  }, error = function(e) {
    message("Error accessing SingleR_labels: ", e$message)
    cell_type_data <- NULL
  })
}

# Third priority: Current Idents
if(is.null(cell_type_data)) {
  tryCatch({
    cell_type_data <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
    message("Using Idents as cell_type")
  }, error = function(e) {
    message("Error accessing Idents: ", e$message)
    # Absolute final fallback
    cell_type_data <- paste0("Cluster_", Peripheral_Blood_Mononuclear_Cells$seurat_clusters)
    message("Using cluster numbers as final fallback")
  })
}

# Step 3: Verify we have valid cell type data
if(is.null(cell_type_data) || length(cell_type_data) != ncol(Peripheral_Blood_Mononuclear_Cells)) {
  message("ERROR: Failed to extract valid cell type data, using emergency fallback")
  message("错误：提取有效细胞类型数据失败，使用紧急回退")
  cell_type_data <- paste0("Cell_", 1:ncol(Peripheral_Blood_Mononuclear_Cells))
}

message("Final cell type data summary:")
message("最终细胞类型数据摘要:")
message("Length: ", length(cell_type_data))
message("Unique types: ", length(unique(cell_type_data)))
message("First 5 types: ", paste(head(cell_type_data, 5), collapse = ", "))

# Create annotation results with robust error handling
annotation_results <- tryCatch({
  data.frame(
    Cell_ID = colnames(Peripheral_Blood_Mononuclear_Cells),
    Cluster = as.character(Peripheral_Blood_Mononuclear_Cells$seurat_clusters),
    Cell_Type = cell_type_data,
    stringsAsFactors = FALSE
  )
}, error = function(e) {
  message("Error creating annotation_results data frame: ", e$message)
  # Emergency data frame creation
  data.frame(
    Cell_ID = colnames(Peripheral_Blood_Mononuclear_Cells),
    Cluster = paste0("Cluster_", 1:ncol(Peripheral_Blood_Mononuclear_Cells)),
    Cell_Type = paste0("Type_", 1:ncol(Peripheral_Blood_Mononuclear_Cells)),
    stringsAsFactors = FALSE
  )
})

write.csv(annotation_results,
          file.path(output_dirs$cell_annotation, "02_Cell_Type_Annotations.csv"),
          row.names = FALSE)

# Add group information
Type <- ifelse(grepl("^control\\.", colnames(Peripheral_Blood_Mononuclear_Cells)), "control", "Disease")
names(Type) <- colnames(Peripheral_Blood_Mononuclear_Cells)
Peripheral_Blood_Mononuclear_Cells <- AddMetaData(
  object = Peripheral_Blood_Mononuclear_Cells,
  metadata = Type,
  col.name = "Type"
)

# Group comparison visualization
pdf(file = file.path(output_dirs$cell_annotation, "03_Group_Comparison.pdf"),
    width = 14, height = 6)
print(
  DimPlot(
    Peripheral_Blood_Mononuclear_Cells,
    reduction = "umap",
    pt.size = 1,
    label = TRUE,
    split.by = "Type"
  ) +
  ggtitle("Cell Types: control vs Disease") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  )
)
dev.off()
# ==============================================================================
# 10. Marker Gene Analysis
# ==============================================================================

message("Performing marker gene analysis")

# Find marker genes for all clusters
markers <- FindAllMarkers(
  object = Peripheral_Blood_Mononuclear_Cells,
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = analysis_params$logFCfilter
)

# Filter significant marker genes
sig_markers <- markers[
  (abs(markers$avg_log2FC) > analysis_params$logFCfilter &
   markers$p_val_adj < analysis_params$adjPvalFilter),
]

# Save results
write.table(
  sig_markers,
  file = file.path(output_dirs$marker_genes, "01_Significant_Marker_Genes.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message(sprintf("Found %d significant marker genes", nrow(sig_markers)))
message(sprintf("找到 %d 个显著marker基因", nrow(sig_markers)))

# Heatmap visualization
if(nrow(sig_markers) > 0) {
  top10 <- sig_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = 10, wt = avg_log2FC)

  pdf(file = file.path(output_dirs$heatmaps, "01_Marker_Genes_Heatmap.pdf"),
      width = 12, height = 9)
  print(
    DoHeatmap(
      object = Peripheral_Blood_Mononuclear_Cells,
      features = top10$gene
    ) +
    NoLegend()
  )
  dev.off()
}

# ==============================================================================
# 11. Differential Expression Analysis
# ==============================================================================

message("Performing differential expression analysis between groups...")
message("进行组间差异表达分析...")

# Get cell type information
cellAnn <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
names(cellAnn) <- colnames(Peripheral_Blood_Mononuclear_Cells)

# Add combined labels
groups <- paste0(Type, "_", cellAnn)
names(groups) <- colnames(Peripheral_Blood_Mononuclear_Cells)
Peripheral_Blood_Mononuclear_Cells <- AddMetaData(
  object = Peripheral_Blood_Mononuclear_Cells,
  metadata = groups,
  col.name = "group"
)

# Differential analysis for each cell type
unique_cell_types <- unique(cellAnn)
differential_results <- list()

for (cellName in unique_cell_types) {
  con_name <- paste0("control_", cellName)
  dis_name <- paste0("Disease_", cellName)

  # Check if labels exist
  if (!(con_name %in% Peripheral_Blood_Mononuclear_Cells$group) |
      !(dis_name %in% Peripheral_Blood_Mononuclear_Cells$group)) {
    message(sprintf("Skipping %s - missing control or disease group", cellName))
    message(sprintf("跳过 %s - 缺少对照组或疾病组", cellName))
    next
  }

  # Get cells - 使用元数据直接筛选，避免 WhichCells expression 参数的兼容性问题
  con_cells <- colnames(Peripheral_Blood_Mononuclear_Cells)[
    Peripheral_Blood_Mononuclear_Cells$group == con_name
  ]
  dis_cells <- colnames(Peripheral_Blood_Mononuclear_Cells)[
    Peripheral_Blood_Mononuclear_Cells$group == dis_name
  ]

  if (length(con_cells) >= 3 & length(dis_cells) >= 3) {
    message(sprintf("Analyzing %s - control: %d cells, Disease: %d cells",
                    cellName, length(con_cells), length(dis_cells)))
    message(sprintf("分析 %s - 对照组: %d 细胞, 疾病组: %d 细胞",
                    cellName, length(con_cells), length(dis_cells)))

    tryCatch({
      pbmc_markers <- FindMarkers(
        Peripheral_Blood_Mononuclear_Cells,
        ident.1 = dis_cells,
        ident.2 = con_cells,
        group.by = 'group',
        logfc.threshold = 0.1
      )

      sig_markers_group <- pbmc_markers[
        (abs(pbmc_markers$avg_log2FC) > analysis_params$logFCfilter &
         pbmc_markers$p_val_adj < analysis_params$adjPvalFilter),
      ]

      if (nrow(sig_markers_group) > 0) {
        sig_markers_group <- cbind(Gene = rownames(sig_markers_group), sig_markers_group)

        # Clean special characters in filename
        safe_cell_name <- gsub("[^A-Za-z0-9_]", "_", cellName)

        write.table(
          sig_markers_group,
          file = file.path(output_dirs$differential_analysis, paste0(safe_cell_name, "_diffGenes.txt")),
          sep = "\t",
          row.names = FALSE,
          quote = FALSE
        )

        differential_results[[cellName]] <- sig_markers_group
        message(sprintf("Saved %d differential genes for %s", nrow(sig_markers_group), cellName))
        message(sprintf("保存 %s 的 %d 个差异基因", cellName, nrow(sig_markers_group)))
      } else {
        message(sprintf("No significant differential genes found for %s", cellName))
        message(sprintf("未找到 %s 的显著差异基因", cellName))
      }
    }, error = function(e) {
      message(sprintf("Error analyzing %s: %s", cellName, e$message))
      message(sprintf("分析 %s 时出错: %s", cellName, e$message))
    })
  } else {
    message(sprintf("Skipping %s - insufficient cells", cellName))
    message(sprintf("跳过 %s - 细胞数量不足", cellName))
  }
}

# ==============================================================================
# 12. Cell Proportion Analysis
# ==============================================================================

message("Performing cell proportion analysis")

# Calculate cell type proportions
cell_type_counts <- table(Idents(Peripheral_Blood_Mononuclear_Cells),
                          Peripheral_Blood_Mononuclear_Cells$Type)

# Calculate proportions
prop_within_group <- prop.table(cell_type_counts, margin = 2) * 100

# Create proportion data frame
proportion_data <- as.data.frame(cell_type_counts)
colnames(proportion_data) <- c("CellType", "Group", "Count")

# Calculate proportions within each group
proportion_summary <- proportion_data %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(
    Total = sum(Count),
    Proportion = Count / Total  # 计算实际比例
  ) %>%
  dplyr::ungroup()

# Save proportion analysis results
write.csv(proportion_summary,
          file.path(output_dirs$proportion_analysis, "01_Cell_Type_Proportions.csv"),
          row.names = FALSE)

# Generate stacked barplot
pdf(file.path(output_dirs$proportion_analysis, "02_Cell_Proportion_Barplot.pdf"),
    width = 8, height = 6)

# Set cell type colors
unique_celltypes <- unique(proportion_summary$CellType)
nColors <- length(unique_celltypes)
celltype_colors <- colorRampPalette(brewer.pal(min(12, nColors), "Set3"))(nColors)
names(celltype_colors) <- unique_celltypes

p <- ggplot(proportion_summary, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", size = 0.3) +
  scale_fill_manual(values = celltype_colors, name = "Cell Type") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(title = "Cell Type Proportions: control vs Disease", y = "Proportion")

print(p)
dev.off()

# ==============================================================================
# 12.2 Sample-wise Cell Type Distribution (Similar to  style)

# ==============================================================================

message("Generating sample-wise cell type distribution plot...")
message("生成按样本的细胞类型分布图...")

# Calculate cell type counts per sample (using orig.ident as sample identifier)
# 计算每个样本中各细胞类型的数量
sample_celltype_counts <- table(
  Idents(Peripheral_Blood_Mononuclear_Cells),
  Peripheral_Blood_Mononuclear_Cells$orig.ident
)

# Convert to data frame / 转换为数据框
sample_proportion_data <- as.data.frame(sample_celltype_counts)
colnames(sample_proportion_data) <- c("CellType", "Sample", "Count")

# Calculate proportions within each sample / 计算每个样本内的比例
sample_proportion_summary <- sample_proportion_data %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(
    Total = sum(Count),
    Proportion = Count / Total,
    Percentage = Proportion * 100
  ) %>%
  dplyr::ungroup()

# Add group information (control vs Disease) based on sample name
# 根据样本名添加分组信息
sample_proportion_summary <- sample_proportion_summary %>%
  dplyr::mutate(
    Group = case_when(
      grepl("control|ctrl|normal|healthy", Sample, ignore.case = TRUE) ~ "Control",
      grepl("disease|treat|tumor|patient", Sample, ignore.case = TRUE) ~ "Disease",
      TRUE ~ Peripheral_Blood_Mononuclear_Cells$Type[match(Sample, Peripheral_Blood_Mononuclear_Cells$orig.ident)]
    )
  )

# If Group is still NA, try to get from metadata
# 如果Group仍然是NA，尝试从元数据获取
if (any(is.na(sample_proportion_summary$Group))) {
  # Create a mapping from sample to group
  sample_group_mapping <- data.frame(
    Sample = unique(Peripheral_Blood_Mononuclear_Cells$orig.ident),
    Group = sapply(unique(Peripheral_Blood_Mononuclear_Cells$orig.ident), function(s) {
      idx <- which(Peripheral_Blood_Mononuclear_Cells$orig.ident == s)[1]
      if (!is.null(Peripheral_Blood_Mononuclear_Cells$Type[idx])) {
        return(as.character(Peripheral_Blood_Mononuclear_Cells$Type[idx]))
      } else {
        return("Unknown")
      }
    })
  )

  sample_proportion_summary <- sample_proportion_summary %>%
    dplyr::select(-Group) %>%
    dplyr::left_join(sample_group_mapping, by = "Sample")
}

# Save sample-wise proportion data / 保存按样本的比例数据
write.csv(sample_proportion_summary,
          file.path(output_dirs$proportion_analysis, "03_Sample_Cell_Type_Proportions.csv"),
          row.names = FALSE)

# Set consistent cell type colors / 设置一致的细胞类型颜色
unique_celltypes_sample <- unique(sample_proportion_summary$CellType)
nColors_sample <- length(unique_celltypes_sample)

# Use a nice color palette / 使用美观的配色方案
if (nColors_sample <= 10) {
  celltype_colors_sample <- brewer.pal(max(3, nColors_sample), "Set3")[1:nColors_sample]
} else {
  celltype_colors_sample <- colorRampPalette(brewer.pal(12, "Set3"))(nColors_sample)
}
names(celltype_colors_sample) <- unique_celltypes_sample

# Order samples by group for better visualization / 按分组排序样本以便更好的可视化
sample_order <- sample_proportion_summary %>%
  dplyr::select(Sample, Group) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Group, Sample)

sample_proportion_summary$Sample <- factor(
  sample_proportion_summary$Sample,
  levels = sample_order$Sample
)

# ============================================
# Plot: Sample-wise Distribution with Group Annotation (Enhanced Version)
# 图：带分组注释的样本分布图（增强版）
# ============================================

pdf(file.path(output_dirs$proportion_analysis, "04_Sample_CellType_Distribution_GroupAnnotated.pdf"),
    width = max(12, length(unique(sample_proportion_summary$Sample)) * 0.7),
    height = 9)

# Create group annotation bar / 创建分组注释条
group_colors <- c("Control" = "#4DBBD5", "Disease" = "#E64B35",
                  "control" = "#4DBBD5", "Disease" = "#E64B35")

# Main stacked bar plot / 主累积条形图
p_main <- ggplot(sample_proportion_summary,
                  aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.85, color = "white", size = 0.15) +
  scale_fill_manual(values = celltype_colors_sample, name = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 8, angle = 60, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(10, 10, 10, 30)
  ) +
  labs(
    title = "Cell Type Distribution Across Samples by Group",
    y = "Cell percent ratio"
  )

print(p_main)
dev.off()

message("Group-annotated sample distribution plot saved!")
message("带分组注释的样本分布图已保存!")

# ============================================
# Plot 3: Cell Type Proportion Summary Statistics
# 图3：细胞类型比例统计汇总
# ============================================

# Calculate summary statistics per cell type and group
# 计算每个细胞类型和分组的统计汇总
celltype_group_summary <- sample_proportion_summary %>%
  dplyr::group_by(CellType, Group) %>%
  dplyr::summarise(
    Mean_Proportion = mean(Proportion),
    SD_Proportion = sd(Proportion),
    Median_Proportion = median(Proportion),
    Min_Proportion = min(Proportion),
    Max_Proportion = max(Proportion),
    N_Samples = n(),
    .groups = "drop"
  )

# Save summary statistics / 保存统计汇总
write.csv(celltype_group_summary,
          file.path(output_dirs$proportion_analysis, "06_CellType_Group_Summary_Statistics.csv"),
          row.names = FALSE)

# Create comparison plot / 创建对比图
pdf(file.path(output_dirs$proportion_analysis, "07_CellType_Proportion_Comparison.pdf"),
    width = 12, height = 8)

p_comparison <- ggplot(celltype_group_summary,
                        aes(x = reorder(CellType, -Mean_Proportion),
                            y = Mean_Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", size = 0.3) +
  geom_errorbar(aes(ymin = Mean_Proportion - SD_Proportion,
                    ymax = Mean_Proportion + SD_Proportion),
                position = position_dodge(width = 0.8), width = 0.25) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Disease" = "#E64B35",
                               "control" = "#4DBBD5", "Disease" = "#E64B35")) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.1))) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    title = "Cell Type Proportion Comparison Between Groups",
    y = "Mean Proportion (± SD)"
  )

print(p_comparison)
dev.off()

message("Cell type proportion comparison plot saved!")
message("细胞类型比例对比图已保存!")

# ==============================================================================
# 13. Gene Expression Visualization for Genes of Interest
# ==============================================================================

message("Generating expression plots for genes of interest...")
message("生成感兴趣基因表达图...")

# Check if genes exist
available_genes <- showGenes[showGenes %in% rownames(GetAssayData(Peripheral_Blood_Mononuclear_Cells, assay = "RNA", layer = "data"))]
missing_genes <- setdiff(showGenes, available_genes)

if (length(missing_genes) > 0) {
  message(paste("The following genes were not found in the dataset:", paste(missing_genes, collapse = ", ")))
  message(paste("以下基因未在数据集中找到:", paste(missing_genes, collapse = ", ")))
}

if (length(available_genes) > 0) {
  # Violin plot
  pdf(file = file.path(output_dirs$gene_expression, "01_ShowGenes_Violin.pdf"),
      width = 15, height = 10)
  print(
    VlnPlot(
      object = Peripheral_Blood_Mononuclear_Cells,
      features = available_genes,
      ncol = min(3, length(available_genes))
    ) +
    ggtitle("Gene Expression Distribution Across Cell Types")
  )
  dev.off()

  # Feature plot
  pdf(file = file.path(output_dirs$gene_expression, "02_ShowGenes_Feature.pdf"),
      width = 15, height = 12)
  print(
    FeaturePlot(
      object = Peripheral_Blood_Mononuclear_Cells,
      features = available_genes,
      cols = c("#E5E5E5", "#FF0000"),
      reduction = "umap",
      ncol = min(3, length(available_genes))
    ) & theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold")
    )
  )
  dev.off()

  # Dot plot - 修复版 (避免细胞类型名称重叠)
  pdf(file = file.path(output_dirs$dotplot_analysis, "01_ShowGenes_DotPlot.pdf"),
      width = 14, height = 9)
  print(
    DotPlot(
      object = Peripheral_Blood_Mononuclear_Cells,
      features = available_genes
    ) +
    coord_flip() +
    theme_minimal(base_size = 11) +
    theme(
      axis.line = element_line(color = "black"),
      # coord_flip() 后: axis.text.x 控制 Cell Types (原始y轴), axis.text.y 控制 Genes (原始x轴)
      axis.text.x = element_text(color = "black", angle = 30, hjust = 1, vjust = 1,
                                  size = 10, face = "bold"),
      axis.text.y = element_text(color = "black", size = 11, face = "italic"),
      axis.title = element_text(color = "black", face = "bold", size = 13),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    labs(title = "Expression of Interest Genes Across Cell Types",
         x = "Genes", y = "Cell Types")
  )
  dev.off()

  # 同步保存 PNG 版本以便查看
  png(file = file.path(output_dirs$dotplot_analysis, "01_ShowGenes_DotPlot.png"),
      width = 14, height = 9, units = "in", res = 300)
  print(
    DotPlot(
      object = Peripheral_Blood_Mononuclear_Cells,
      features = available_genes
    ) +
    coord_flip() +
    theme_minimal(base_size = 11) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(color = "black", angle = 30, hjust = 1, vjust = 1,
                                  size = 10, face = "bold"),
      axis.text.y = element_text(color = "black", size = 11, face = "italic"),
      axis.title = element_text(color = "black", face = "bold", size = 13),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    labs(title = "Expression of Interest Genes Across Cell Types",
         x = "Genes", y = "Cell Types")
  )
  dev.off()
}

# ==============================================================================
# 13. Advanced BoxPlot Analysis for Show Genes
# ==============================================================================

message("Performing advanced boxplot analysis for genes of interest...")
message("对感兴趣基因进行高级箱线图分析...")

# Check if genes exist
available_genes <- showGenes[showGenes %in% rownames(GetAssayData(Peripheral_Blood_Mononuclear_Cells, assay = "RNA", layer = "data"))]
missing_genes <- setdiff(showGenes, available_genes)

if (length(missing_genes) > 0) {
  message(paste("The following genes were not found in the dataset:", paste(missing_genes, collapse = ", ")))
  message(paste("以下基因未在数据集中找到:", paste(missing_genes, collapse = ", ")))
}

if (length(available_genes) > 0) {
  message(sprintf("Starting analysis of %d available genes", length(available_genes)))
  message(sprintf("开始分析 %d 个可用基因", length(available_genes)))

  # Advanced BoxPlot Analysis Function
  create_advanced_boxplots <- function(seurat_obj, genes, cell_ann, type_info) {

    # Get expression data
    expr_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

    # Create data frame
    plot_data <- data.frame(
      Cell = colnames(seurat_obj),
      CellType = cell_ann,
      Group = type_info,
      stringsAsFactors = FALSE
    )

    # Add gene expression data
    for (gene in genes) {
      plot_data[[gene]] <- as.numeric(expr_data[gene, ])
    }

    # Statistical results storage
    all_stats <- data.frame()

    # Generate boxplots for each gene
    for (gene in genes) {
      message(sprintf("Processing gene: %s", gene))
      message(sprintf("正在处理基因: %s", gene))

      # Prepare current gene data
      current_data <- plot_data[, c("Cell", "CellType", "Group", gene)]
      colnames(current_data)[4] <- "Expression"

      # Calculate expression range for y-axis
      expr_values <- current_data$Expression
      y_max <- quantile(expr_values, 0.95, na.rm = TRUE)
      y_min <- 0

      # Calculate statistics for each cell type
      cell_types <- unique(current_data$CellType)
      gene_stats <- data.frame()

      for (ct in cell_types) {
        ct_data <- current_data[current_data$CellType == ct, ]

        # Check if sufficient cells for statistics
        control_cells <- sum(ct_data$Group == "control")
        disease_cells <- sum(ct_data$Group == "Disease")

        if (control_cells >= 3 && disease_cells >= 3) {
          # Perform Wilcoxon rank-sum test
          tryCatch({
            test_result <- wilcox.test(
              Expression ~ Group,
              data = ct_data,
              alternative = "two.sided"
            )

            p_value <- test_result$p.value

            # Calculate means and medians
            control_mean <- mean(ct_data$Expression[ct_data$Group == "control"], na.rm = TRUE)
            disease_mean <- mean(ct_data$Expression[ct_data$Group == "Disease"], na.rm = TRUE)
            control_median <- median(ct_data$Expression[ct_data$Group == "control"], na.rm = TRUE)
            disease_median <- median(ct_data$Expression[ct_data$Group == "Disease"], na.rm = TRUE)

            # Calculate fold change
            fold_change <- ifelse(control_mean > 0, disease_mean / control_mean, NA)
            log2_fc <- ifelse(control_mean > 0, log2(disease_mean / control_mean), NA)

            # Add to statistical results
            stat_row <- data.frame(
              Gene = gene,
              CellType = ct,
              control_Mean = control_mean,
              Disease_Mean = disease_mean,
              control_Median = control_median,
              Disease_Median = disease_median,
              Fold_Change = fold_change,
              Log2_FC = log2_fc,
              P_value = p_value,
              P_adj = p.adjust(p_value, method = "BH"),
              Significant = p_value < 0.05,
              P_symbol = case_when(
                p_value < 0.001 ~ "***",
                p_value < 0.01 ~ "**",
                p_value < 0.05 ~ "*",
                TRUE ~ ""
              ),
              control_Cells = control_cells,
              Disease_Cells = disease_cells
            )

            gene_stats <- rbind(gene_stats, stat_row)

          }, error = function(e) {
            message(sprintf("Statistical test failed - Gene: %s, Cell Type: %s, Error: %s",
                            gene, ct, e$message))
            message(sprintf("统计检验失败 - 基因: %s, 细胞类型: %s, 错误: %s",
                            gene, ct, e$message))
          })
        }
      }

      # Save individual gene statistics
      if (nrow(gene_stats) > 0) {
        write.csv(gene_stats,
                  file.path(output_dirs$boxplot_analysis, "ShowGenes_Statistics", paste0(gene, "_boxplot_statistics.csv")),
                  row.names = FALSE)
      }

      # Create individual boxplot for current gene
      tryCatch({
        # Filter cell types with sufficient cells
        valid_celltypes <- gene_stats$CellType[gene_stats$control_Cells >= 3 & gene_stats$Disease_Cells >= 3]
        plot_data_filtered <- current_data[current_data$CellType %in% valid_celltypes, ]

        if (nrow(plot_data_filtered) > 0) {
          p <- ggplot(plot_data_filtered, aes(x = CellType, y = Expression, fill = Group)) +
            geom_boxplot(alpha = 0.8,
                         position = position_dodge(width = 0.8),
                         outlier.size = 0.5) +
            geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
                       alpha = 0.4, size = 0.3) +
            scale_fill_manual(values = c("control" = "#4CAF50", "Disease" = "#F44336"),
                             name = "Group") +
            labs(title = paste0(gene, " Expression by Cell Type and Group"),
                 subtitle = paste("n =", nrow(plot_data_filtered), "cells"),
                 x = "Cell Type",
                 y = "Expression Level") +
            theme_minimal(base_size = 12) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
              axis.text.y = element_text(face = "bold"),
              axis.title = element_text(face = "bold", size = 14),
              plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 12),
              legend.position = "top",
              legend.title = element_text(face = "bold"),
              legend.text = element_text(face = "bold"),
              panel.background = element_rect(fill = "white"),
              plot.background = element_rect(fill = "white")
            ) +
            ylim(y_min, y_max * 1.1) +
            NoGrid()

          # Add significance annotations
          if (nrow(gene_stats) > 0) {
            for (i in 1:nrow(gene_stats)) {
              if (!is.na(gene_stats$Significant[i]) && gene_stats$Significant[i]) {
                ct <- gene_stats$CellType[i]
                p_symbol <- gene_stats$P_symbol[i]
                y_pos <- y_max * 1.05

                p <- p + annotate("text",
                                  x = ct,
                                  y = y_pos,
                                  label = p_symbol,
                                  size = 5,
                                  fontface = "bold",
                                  color = "red")
              }
            }
          }

          # Save individual plots
          pdf_file <- file.path(output_dirs$boxplot_analysis, "ShowGenes_Individual_Plots", paste0(gene, "_BoxPlot_controlVsDisease.pdf"))
          ggsave(pdf_file, plot = p, width = 12, height = 8, device = "pdf")

          png_file <- file.path(output_dirs$boxplot_analysis, "ShowGenes_Individual_Plots", paste0(gene, "_BoxPlot_controlVsDisease.png"))
          ggsave(png_file, plot = p, width = 12, height = 8, device = "png", dpi = 300)

          message(sprintf("Individual boxplot saved for gene: %s", gene))
          message(sprintf("基因 %s 的单独箱线图已保存", gene))
        }

      }, error = function(e) {
        message(sprintf("Error creating boxplot for gene %s: %s", gene, e$message))
        message(sprintf("为基因 %s 创建箱线图时出错: %s", gene, e$message))
      })

      # Add current gene stats to overall stats
      all_stats <- rbind(all_stats, gene_stats)
    }

    # Save overall statistics
    write.csv(all_stats, file.path(output_dirs$boxplot_analysis, "ShowGenes_Statistics", "All_Genes_BoxPlot_Statistics.csv"), row.names = FALSE)

    return(all_stats)
  }

  # Run advanced boxplot analysis
  message("Running advanced boxplot analysis...")
  message("运行高级箱线图分析...")

  # Get cell annotations
  cell_annotations <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
  names(cell_annotations) <- colnames(Peripheral_Blood_Mononuclear_Cells)

  # Get group information
  group_info <- Peripheral_Blood_Mononuclear_Cells$Type
  names(group_info) <- colnames(Peripheral_Blood_Mononuclear_Cells)

  # Run the analysis
  boxplot_stats <- create_advanced_boxplots(
    Peripheral_Blood_Mononuclear_Cells,
    available_genes,
    cell_annotations,
    group_info
  )

  # Create comprehensive boxplot
  message("Creating comprehensive boxplot...")
  message("创建综合箱线图...")

  if (nrow(boxplot_stats) > 0) {
    # Create comprehensive plot data
    expr_data <- GetAssayData(Peripheral_Blood_Mononuclear_Cells, assay = "RNA", layer = "data")

    comp_data <- data.frame()
    for (gene in available_genes) {
      gene_data <- data.frame(
        Gene = gene,
        CellType = cell_annotations,
        Group = group_info,
        Expression = as.numeric(expr_data[gene, ]),
        stringsAsFactors = FALSE
      )
      comp_data <- rbind(comp_data, gene_data)
    }

    # Create comprehensive boxplot
    comp_plot <- ggplot(comp_data, aes(x = CellType, y = Expression, fill = Group)) +
      geom_boxplot(alpha = 0.8,
                   position = position_dodge(width = 0.8),
                   outlier.size = 0.3) +
      scale_fill_manual(values = c("control" = "#4CAF50", "Disease" = "#F44336"),
                       name = "Group") +
      facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
      labs(title = "Expression of Interest Genes: control vs Disease",
           subtitle = paste("Genes:", paste(available_genes, collapse = ", ")),
           x = "Cell Type",
           y = "Expression Level") +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 10),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")
      ) +
      NoGrid()

    # Save comprehensive plots
    pdf_file <- file.path(output_dirs$boxplot_analysis, "ShowGenes_BoxPlots", "Comprehensive_ShowGenes_BoxPlots.pdf")
    ggsave(pdf_file, plot = comp_plot, width = 16, height = 12, device = "pdf")

    png_file <- file.path(output_dirs$boxplot_analysis, "ShowGenes_BoxPlots", "Comprehensive_ShowGenes_BoxPlots.png")
    ggsave(png_file, plot = comp_plot, width = 16, height = 12, device = "png", dpi = 300)

    message("Comprehensive boxplot saved")
    message("综合箱线图已保存")

    # Generate summary report
    gene_summary <- boxplot_stats %>%
      group_by(Gene) %>%
      summarise(
        CellTypes_Tested = n(),
        Significant_CellTypes = sum(Significant, na.rm = TRUE),
        Max_Log2FC = max(abs(Log2_FC), na.rm = TRUE),
        Min_Pvalue = min(P_value, na.rm = TRUE),
        .groups = "drop"
      )

    celltype_summary <- boxplot_stats %>%
      group_by(CellType) %>%
      summarise(
        Genes_Tested = n(),
        Significant_Genes = sum(Significant, na.rm = TRUE),
        Avg_Log2FC = mean(Log2_FC, na.rm = TRUE),
        .groups = "drop"
      )

    significant_combinations <- boxplot_stats %>%
      dplyr::filter(Significant) %>%
      arrange(P_value) %>%
      dplyr::select(Gene, CellType, Log2_FC, P_value, P_symbol)

    # Save summary reports
    write.csv(gene_summary, file.path(output_dirs$boxplot_analysis, "ShowGenes_Statistics", "Gene_Summary_BoxPlot.csv"), row.names = FALSE)
    write.csv(celltype_summary, file.path(output_dirs$boxplot_analysis, "ShowGenes_Statistics", "CellType_Summary_BoxPlot.csv"), row.names = FALSE)
    write.csv(significant_combinations, file.path(output_dirs$boxplot_analysis, "ShowGenes_Statistics", "Significant_Combinations_BoxPlot.csv"), row.names = FALSE)

    message("BoxPlot analysis completed successfully!")
    message("箱线图分析成功完成！")
    message(sprintf("Total significant gene-celltype combinations: %d", nrow(significant_combinations)))
    message(sprintf("显著的基因-细胞类型组合总数: %d", nrow(significant_combinations)))
  }

} else {
  message("No genes available for boxplot analysis")
  message("没有可用于箱线图分析的基因")
}

# ==============================================================================
# 14. Individual Gene UMAP Distribution Analysis
# ==============================================================================

message("Performing individual gene UMAP distribution analysis...")
message("进行单基因UMAP分布分析...")

# Create output directory for individual gene analysis
individual_gene_dir <- "17.Individual_Gene_UMAP"
output_dirs$individual_gene_umap <- individual_gene_dir

# Create main directory and subdirectories
full_path <- file.path(workDir, individual_gene_dir)
if (!dir.exists(full_path)) {
  dir.create(full_path, recursive = TRUE)
  message(sprintf("Created directory: %s", individual_gene_dir))
}

# Create subdirectories for each gene
gene_subdirs <- c("Individual_FeaturePlots", "Combined_Analysis", "Gene_Statistics")
for (subdir in gene_subdirs) {
  subdir_path <- file.path(full_path, subdir)
  if (!dir.exists(subdir_path)) {
    dir.create(subdir_path, recursive = TRUE)
  }
}

# Check if genes exist
available_genes <- showGenes[showGenes %in% rownames(GetAssayData(Peripheral_Blood_Mononuclear_Cells, assay = "RNA", layer = "data"))]
missing_genes <- setdiff(showGenes, available_genes)

if (length(missing_genes) > 0) {
  message(paste("The following genes were not found in the dataset:", paste(missing_genes, collapse = ", ")))
  message(paste("以下基因未在数据集中找到:", paste(missing_genes, collapse = ", ")))
}

if (length(available_genes) > 0) {
  message(sprintf("Starting individual gene UMAP analysis for %d available genes", length(available_genes)))
  message(sprintf("开始为 %d 个可用基因进行单基因UMAP分析", length(available_genes)))

  # Function to create individual gene UMAP plots
  create_individual_gene_plots <- function(seurat_obj, gene, output_dir) {

    # 1. Basic FeaturePlot
    message(sprintf("Creating basic FeaturePlot for gene: %s", gene))

    pdf(file.path(output_dir, "Individual_FeaturePlots", paste0(gene, "_UMAP_FeaturePlot.pdf")),
        width = 8, height = 6)
    print(
      FeaturePlot(
        seurat_obj,
        features = gene,
        reduction = "umap",
        cols = c("lightgrey", "#FF0000"),
        pt.size = 1.2,
        label = TRUE,
        label.size = 3
      ) +
      ggtitle(paste(gene, "Expression (UMAP)")) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
    dev.off()

    # 2. Combined analysis with cell type annotation (UMAP + Legend style)
    message(sprintf("Creating combined analysis plot for gene: %s", gene))

    # Get UMAP embeddings
    umap_df <- as.data.frame(Embeddings(seurat_obj, "umap"))
    colnames(umap_df) <- c("UMAP_1", "UMAP_2")

    # Add metadata
    umap_df$Cluster <- as.character(seurat_obj$seurat_clusters)
    umap_df$CellType <- as.character(Idents(seurat_obj))
    umap_df$Group <- seurat_obj$Type
    umap_df$Expression <- FetchData(seurat_obj, vars = gene)[,1]

    # Calculate cluster centers for labeling
    label_df <- umap_df %>%
      group_by(Cluster, CellType) %>%
      summarise(
        UMAP_1 = mean(UMAP_1),
        UMAP_2 = mean(UMAP_2),
        .groups = "drop"
      )

    # Calculate cluster statistics for bar plot
    cluster_stats <- umap_df %>%
      group_by(Cluster, CellType) %>%
      summarise(
        AvgExpr = mean(Expression),
        .groups = "drop"
      ) %>%
      arrange(as.numeric(Cluster))

    # Create color scheme for clusters
    cluster_levels <- sort(unique(umap_df$Cluster))
    cluster_colors <- setNames(colorRampPalette(brewer.pal(min(12, length(cluster_levels)), "Set3"))(length(cluster_levels)), cluster_levels)

    # UMAP main plot (similar to reference style)
    p_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = factor(Cluster)), size = 1.2, alpha = 0.7) +
      scale_color_manual(values = cluster_colors, name = "Cluster") +
      geom_text(
        data = label_df,
        aes(label = paste0("C", Cluster, "\\n", CellType)),
        color = "black",
        size = 3.5,
        fontface = "bold"
      ) +
      theme_classic(base_size = 14) +
      ggtitle(gene) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", face = "bold")
      )

    # Bar plot (legend) showing average expression
    cluster_stats$ClusterLabel <- paste0("C", cluster_stats$Cluster)
    cluster_stats$ClusterLabel <- factor(cluster_stats$ClusterLabel, levels = paste0("C", cluster_levels))
    max_expr <- max(cluster_stats$AvgExpr)
    offset <- max_expr * 0.05

    p_legend <- ggplot(cluster_stats, aes(x = ClusterLabel, y = AvgExpr, fill = Cluster)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = cluster_colors, guide = "none") +
      geom_text(
        aes(x = ClusterLabel, y = AvgExpr + offset, label = CellType),
        hjust = 0, vjust = 0.5, size = 4, fontface = "plain", color = "black"
      ) +
      coord_flip(clip = "off") +
      labs(x = NULL, y = "Average Expression", title = "") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12, face = "bold"),
        plot.margin = ggplot2::margin(5, 80, 5, 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
      )

    # Combined plot using patchwork (like reference style)
    combined_plot <- p_umap + p_legend + plot_layout(widths = c(2.5, 1.2))

    pdf(file.path(output_dir, "Combined_Analysis", paste0(gene, "_UMAP_Combined_Analysis.pdf")),
        width = 12, height = 6)
    print(combined_plot)
    dev.off()

    # Also create individual group comparison plots
    message(sprintf("Creating group comparison plots for gene: %s", gene))

    for (group_name in unique(umap_df$Group)) {
      group_data <- umap_df[umap_df$Group == group_name, ]

      # Calculate group-specific cluster centers
      group_label_df <- group_data %>%
        group_by(Cluster, CellType) %>%
        summarise(
          UMAP_1 = mean(UMAP_1),
          UMAP_2 = mean(UMAP_2),
          .groups = "drop"
        )

      # Group-specific cluster statistics
      group_cluster_stats <- group_data %>%
        group_by(Cluster, CellType) %>%
        summarise(
          AvgExpr = mean(Expression),
          .groups = "drop"
        ) %>%
        arrange(as.numeric(Cluster))

      # UMAP plot for this group
      p_group_umap <- ggplot(group_data, aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(aes(color = factor(Cluster)), size = 1.2, alpha = 0.7) +
        scale_color_manual(values = cluster_colors, name = "Cluster") +
        geom_text(
          data = group_label_df,
          aes(label = paste0("C", Cluster)),
          color = "black",
          size = 3.5,
          fontface = "bold"
        ) +
        theme_classic(base_size = 14) +
        ggtitle(paste(gene, "-", group_name)) +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "none",
          axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black", face = "bold")
        )

      # Bar plot for this group
      if (nrow(group_cluster_stats) > 0) {
        group_cluster_stats$ClusterLabel <- paste0("C", group_cluster_stats$Cluster)
        group_cluster_stats$ClusterLabel <- factor(group_cluster_stats$ClusterLabel,
                                                   levels = paste0("C", cluster_levels))
        group_max_expr <- max(group_cluster_stats$AvgExpr)
        group_offset <- group_max_expr * 0.05

        p_group_legend <- ggplot(group_cluster_stats, aes(x = ClusterLabel, y = AvgExpr, fill = Cluster)) +
          geom_bar(stat = "identity", width = 0.7) +
          scale_fill_manual(values = cluster_colors, guide = "none") +
          geom_text(
            aes(x = ClusterLabel, y = AvgExpr + group_offset, label = CellType),
            hjust = 0, vjust = 0.5, size = 4, fontface = "plain", color = "black"
          ) +
          coord_flip(clip = "off") +
          labs(x = NULL, y = "Average Expression", title = "") +
          theme_minimal(base_size = 12) +
          theme(
            axis.text.y = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 10),
            axis.title.x = element_text(size = 12, face = "bold"),
            plot.margin = ggplot2::margin(5, 80, 5, 5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black")
          )

        # Combined group plot
        group_combined_plot <- p_group_umap + p_group_legend + plot_layout(widths = c(2.5, 1.2))

        pdf(file.path(output_dir, "Combined_Analysis", paste0(gene, "_UMAP_", group_name, "_Analysis.pdf")),
            width = 12, height = 6)
        print(group_combined_plot)
        dev.off()
      }
    }

    # 4. Export expression statistics
    message(sprintf("Exporting statistics for gene: %s", gene))

    # Average expression by cell type and group
    avg_by_celltype_group <- umap_df %>%
      group_by(CellType, Group) %>%
      summarise(
        Cell_Count = n(),
        Mean_Expression = mean(Expression),
        Median_Expression = median(Expression),
        SD_Expression = sd(Expression),
        .groups = "drop"
      ) %>%
      arrange(CellType, Group)

    # Average expression by cluster and group
    avg_by_cluster_group <- umap_df %>%
      group_by(Cluster, Group) %>%
      summarise(
        Cell_Count = n(),
        Mean_Expression = mean(Expression),
        Median_Expression = median(Expression),
        SD_Expression = sd(Expression),
        .groups = "drop"
      ) %>%
      arrange(Cluster, Group)

    # Save statistics
    write.csv(avg_by_celltype_group,
              file.path(output_dir, "Gene_Statistics", paste0(gene, "_Expression_by_CellType_Group.csv")),
              row.names = FALSE)
    write.csv(avg_by_cluster_group,
              file.path(output_dir, "Gene_Statistics", paste0(gene, "_Expression_by_Cluster_Group.csv")),
              row.names = FALSE)

    # Save individual cell expression data
    cell_expr_data <- data.frame(
      Cell_ID = rownames(umap_df),
      Cluster = umap_df$Cluster,
      CellType = umap_df$CellType,
      Group = umap_df$Group,
      UMAP_1 = umap_df$UMAP_1,
      UMAP_2 = umap_df$UMAP_2,
      Expression = umap_df$Expression
    )

    write.csv(cell_expr_data,
              file.path(output_dir, "Gene_Statistics", paste0(gene, "_Individual_Cell_Expression.csv")),
              row.names = FALSE)

    message(sprintf("Completed analysis for gene: %s", gene))
    message(sprintf("基因 %s 分析完成", gene))

    return(list(
      expression_data = umap_df,
      celltype_stats = avg_by_celltype_group,
      cluster_stats = avg_by_cluster_group
    ))
  }

  # Process each gene
  gene_results <- list()
  for (gene in available_genes) {
    tryCatch({
      gene_results[[gene]] <- create_individual_gene_plots(
        Peripheral_Blood_Mononuclear_Cells,
        gene,
        full_path
      )
    }, error = function(e) {
      message(sprintf("Error processing gene %s: %s", gene, e$message))
      message(sprintf("处理基因 %s 时出错: %s", gene, e$message))
    })
  }

  # Create comprehensive summary report
  message("Creating comprehensive gene expression summary...")
  message("创建综合基因表达汇总...")

  if (length(gene_results) > 0) {
    # Combine all cell type statistics
    all_celltype_stats <- data.frame()
    all_cluster_stats <- data.frame()

    for (gene in names(gene_results)) {
      if (!is.null(gene_results[[gene]]$celltype_stats)) {
        celltype_stats <- gene_results[[gene]]$celltype_stats
        celltype_stats$Gene <- gene
        all_celltype_stats <- rbind(all_celltype_stats, celltype_stats)
      }

      if (!is.null(gene_results[[gene]]$cluster_stats)) {
        cluster_stats <- gene_results[[gene]]$cluster_stats
        cluster_stats$Gene <- gene
        all_cluster_stats <- rbind(all_cluster_stats, cluster_stats)
      }
    }

    # Save comprehensive summaries
    write.csv(all_celltype_stats,
              file.path(full_path, "Gene_Statistics", "All_Genes_CellType_Expression_Summary.csv"),
              row.names = FALSE)
    write.csv(all_cluster_stats,
              file.path(full_path, "Gene_Statistics", "All_Genes_Cluster_Expression_Summary.csv"),
              row.names = FALSE)

    # Create summary statistics
    gene_summary <- all_celltype_stats %>%
      group_by(Gene, Group) %>%
      summarise(
        Total_CellTypes = n(),
        Avg_Expression_Across_CellTypes = mean(Mean_Expression, na.rm = TRUE),
        Max_Expression = max(Mean_Expression, na.rm = TRUE),
        Min_Expression = min(Mean_Expression, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = Group,
        values_from = c(Avg_Expression_Across_CellTypes, Max_Expression, Min_Expression),
        names_sep = "_"
      )

    write.csv(gene_summary,
              file.path(full_path, "Gene_Statistics", "Gene_Expression_Summary_by_Group.csv"),
              row.names = FALSE)

    message("Individual gene UMAP distribution analysis completed successfully!")
    message("单基因UMAP分布分析成功完成！")
    message(sprintf("Analyzed %d genes with complete statistics and visualizations", length(available_genes)))
    message(sprintf("分析了 %d 个基因，包含完整的统计和可视化结果", length(available_genes)))
  }

} else {
  message("No genes available for individual UMAP analysis")
  message("没有可用于单基因UMAP分析的基因")
}

# ==============================================================================
# 15. Circular UMAP Visualization (using plot1cell)
# ==============================================================================

message("Generating Circular UMAP plots...")
message("生成环形UMAP图...")

# Check if plot1cell package is available
if (!requireNamespace("plot1cell", quietly = TRUE)) {
  message("Warning: plot1cell package not installed, skipping circular UMAP plots")
  message("警告：plot1cell包未安装，跳过环形UMAP图")
} else {
  tryCatch({
    # =============================================
    # 1. All cells circular UMAP
    # =============================================
    message("Creating circular UMAP for all cells...")
    message("创建所有细胞的环形UMAP图...")

    # Set Idents to cell type
    Idents(Peripheral_Blood_Mononuclear_Cells) <- Peripheral_Blood_Mononuclear_Cells@meta.data$cell_type
    celltype_order_all <- sort(unique(Peripheral_Blood_Mononuclear_Cells@meta.data$cell_type))
    levels(Peripheral_Blood_Mononuclear_Cells) <- celltype_order_all

    # Prepare circular plot data
    circ_data_all <- prepare_circlize_data(Peripheral_Blood_Mononuclear_Cells, scale = 0.8)

    # Generate colors
    set.seed(1234)
    cluster_colors_all <- rand_color(length(levels(Peripheral_Blood_Mononuclear_Cells)))

    # Save PDF
    pdf(file.path(output_dirs$circular_umap, "01_Circular_UMAP_All_Cells.pdf"),
        width = 6, height = 6)

    plot_circlize(circ_data_all,
                  do.label = TRUE,
                  pt.size = 0.1,
                  col.use = cluster_colors_all,
                  bg.color = 'white',
                  kde2d.n = 200,
                  repel = TRUE,
                  label.cex = 0.6)

    dev.off()

    message("All cells circular UMAP saved")
    message("所有细胞环形UMAP图已保存")

    # =============================================
    # 2. control group circular UMAP
    # =============================================
    message("Creating circular UMAP for control group...")
    message("创建对照组环形UMAP图...")

    # Create control subset
    scObject_control <- subset(Peripheral_Blood_Mononuclear_Cells, subset = Type == "control")

    if (ncol(scObject_control) > 0) {
      # Check cell type counts
      celltype_counts <- table(scObject_control@meta.data$cell_type)
      min_cells_per_type <- 5
      valid_celltypes <- names(celltype_counts[celltype_counts >= min_cells_per_type])

      if (length(valid_celltypes) > 0) {
        # Filter to keep valid cell types
        cells_to_keep <- scObject_control@meta.data$cell_type %in% valid_celltypes
        scObject_control_filtered <- scObject_control[, cells_to_keep]

        # Set Idents
        Idents(scObject_control_filtered) <- scObject_control_filtered@meta.data$cell_type
        celltype_order_control <- sort(unique(scObject_control_filtered@meta.data$cell_type))
        levels(scObject_control_filtered) <- celltype_order_control

        # Prepare data
        circ_data_control <- prepare_circlize_data(scObject_control_filtered, scale = 0.8)

        # Generate colors
        set.seed(1234)
        cluster_colors_control <- rand_color(length(levels(scObject_control_filtered)))

        # Save PDF
        pdf(file.path(output_dirs$circular_umap, "02_Circular_UMAP_control.pdf"),
            width = 6, height = 6)

        plot_circlize(circ_data_control,
                      do.label = TRUE,
                      pt.size = 0.1,
                      col.use = cluster_colors_control,
                      bg.color = 'white',
                      kde2d.n = 200,
                      repel = TRUE,
                      label.cex = 0.6)

        dev.off()

        message("control group circular UMAP saved")
        message("对照组环形UMAP图已保存")
      } else {
        message("Warning: Not enough cells per type in control group")
        message("警告：对照组各细胞类型细胞数不足")
      }
    } else {
      message("Warning: No control group cells found")
      message("警告：未找到对照组细胞")
    }

    # =============================================
    # 3. Disease group circular UMAP
    # =============================================
    message("Creating circular UMAP for Disease group...")
    message("创建疾病组环形UMAP图...")

    # Create Disease subset
    scObject_disease <- subset(Peripheral_Blood_Mononuclear_Cells, subset = Type == "Disease")

    if (ncol(scObject_disease) > 0) {
      # Check cell type counts
      celltype_counts <- table(scObject_disease@meta.data$cell_type)
      min_cells_per_type <- 5
      valid_celltypes <- names(celltype_counts[celltype_counts >= min_cells_per_type])

      if (length(valid_celltypes) > 0) {
        # Filter to keep valid cell types
        cells_to_keep <- scObject_disease@meta.data$cell_type %in% valid_celltypes
        scObject_disease_filtered <- scObject_disease[, cells_to_keep]

        # Set Idents
        Idents(scObject_disease_filtered) <- scObject_disease_filtered@meta.data$cell_type
        celltype_order_disease <- sort(unique(scObject_disease_filtered@meta.data$cell_type))
        levels(scObject_disease_filtered) <- celltype_order_disease

        # Prepare data
        circ_data_disease <- prepare_circlize_data(scObject_disease_filtered, scale = 0.8)

        # Generate colors
        set.seed(1234)
        cluster_colors_disease <- rand_color(length(levels(scObject_disease_filtered)))

        # Save PDF
        pdf(file.path(output_dirs$circular_umap, "03_Circular_UMAP_Disease.pdf"),
            width = 6, height = 6)

        plot_circlize(circ_data_disease,
                      do.label = TRUE,
                      pt.size = 0.1,
                      col.use = cluster_colors_disease,
                      bg.color = 'white',
                      kde2d.n = 200,
                      repel = TRUE,
                      label.cex = 0.6)

        dev.off()

        message("Disease group circular UMAP saved")
        message("疾病组环形UMAP图已保存")
      } else {
        message("Warning: Not enough cells per type in Disease group")
        message("警告：疾病组各细胞类型细胞数不足")
      }
    } else {
      message("Warning: No Disease group cells found")
      message("警告：未找到疾病组细胞")
    }

    message("Circular UMAP visualization completed!")
    message("环形UMAP图可视化完成！")

  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message(sprintf("Error in circular UMAP generation: %s", e$message))
    message(sprintf("环形UMAP图生成错误: %s", e$message))
  })
}

# ==============================================================================
# 16. scRNAtoolVis Visualization / scRNAtoolVis可视化
# ==============================================================================

message("Generating scRNAtoolVis visualizations...")
message("生成scRNAtoolVis可视化图形...")

# --- 16.1 Cell Ratio Barplot / 细胞比例柱状图 ---
message("Creating cell ratio barplot...")
message("创建细胞比例柱状图...")

# Add Sample metadata (use orig.ident as sample identifier)
# 添加Sample元数据（使用orig.ident作为样本标识）
if (!"Sample" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data)) {
  Peripheral_Blood_Mononuclear_Cells$Sample <- Peripheral_Blood_Mononuclear_Cells$Type
}

# Ensure cellType column exists for scRNAtoolVis functions
# 确保cellType列存在以供scRNAtoolVis函数使用
if (!"cellType" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data)) {
  Peripheral_Blood_Mononuclear_Cells$cellType <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
}

pdf(file = file.path(output_dirs$scRNAtoolVis, "01_cellRatioPlot.pdf"), width = 6, height = 6)
cellRatioPlot(
  object = Peripheral_Blood_Mononuclear_Cells,
  sample.name = "Sample",
  celltype.name = "cellType",
  col.width = 0.7
)
dev.off()

message("Cell ratio barplot saved / 细胞比例柱状图已保存")

# --- 16.2 Marker Volcano Plot / Marker基因火山图 ---
message("Creating marker volcano plot...")
message("创建Marker基因火山图...")

# Use the markers from FindAllMarkers (section 10)
# 使用FindAllMarkers的结果（第10节）
if (exists("markers") && nrow(markers) > 0) {
  # Determine number of cell types for color palette
  # 根据细胞类型数量生成足够的颜色
  n_celltypes <- length(unique(markers$cluster))
  if (n_celltypes <= 10) {
    volcano_colors <- ggsci::pal_npg()(n_celltypes)
  } else {
    volcano_colors <- colorRampPalette(ggsci::pal_npg()(10))(n_celltypes)
  }

  tryCatch({
    pdf(file = file.path(output_dirs$scRNAtoolVis, "02_markerVolcano.pdf"), width = 12, height = 6)
    p_markerVolcano <- markerVolcano(markers = markers, topn = 5,
                  labelCol = volcano_colors,
                  log2FC = analysis_params$logFCfilter)
    print(p_markerVolcano)
    dev.off()
    message("Marker volcano plot saved / Marker基因火山图已保存")
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message(sprintf("markerVolcano error: %s", e$message))
    message("尝试不指定labelCol重新绑制...")
    tryCatch({
      pdf(file = file.path(output_dirs$scRNAtoolVis, "02_markerVolcano.pdf"), width = 10, height = 6)
      p_markerVolcano2 <- markerVolcano(markers = markers, topn = 5,
                    log2FC = analysis_params$logFCfilter)
      print(p_markerVolcano2)
      dev.off()
      message("Marker volcano plot saved (default colors) / Marker基因火山图已保存（默认颜色）")
    }, error = function(e2) {
      if (dev.cur() > 1) dev.off()
      message(sprintf("markerVolcano retry error: %s", e2$message))
    })
  })

  # Generate enough colors for all cell type clusters in jjVolcano
  # jjVolcano的颜色参数是tile.col，默认只有9色，需要根据cluster数量扩展
  jj_n_celltypes <- length(unique(markers$cluster))
  jj_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(jj_n_celltypes)

  # --- 16.3 jjVolcano Plot 1 (Normal) / jjVolcano火山图1（标准） ---
  message("Creating jjVolcano plot 1 (normal)...")
  message("创建jjVolcano火山图1（标准）...")

  pdf(file = file.path(output_dirs$scRNAtoolVis, "03_jjVolcano_normal.pdf"), width = 12, height = 6)
  print(
    jjVolcano(diffData = markers,
              legend.position = c(0.9, 0.95),
              tile.col = jj_colors,
              log2FC.cutoff = analysis_params$logFCfilter) +
      labs(x = "Cell type")
  )
  dev.off()

  message("jjVolcano plot 1 saved / jjVolcano火山图1已保存")

  # --- 16.4 jjVolcano Plot 2 (Flipped) / jjVolcano火山图2（转置） ---
  message("Creating jjVolcano plot 2 (flipped)...")
  message("创建jjVolcano火山图2（转置）...")

  pdf(file = file.path(output_dirs$scRNAtoolVis, "04_jjVolcano_flipped.pdf"), width = 9, height = 7)
  print(
    jjVolcano(diffData = markers,
              legend.position = c(0.95, 0.9),
              flip = TRUE,
              tile.col = jj_colors,
              log2FC.cutoff = analysis_params$logFCfilter) +
      labs(x = "Cell type")
  )
  dev.off()

  message("jjVolcano plot 2 saved / jjVolcano火山图2已保存")

  # --- 16.5 jjVolcano Plot 3 (Polar/Circular) / jjVolcano火山图3（环形） ---
  message("Creating jjVolcano plot 3 (polar)...")
  message("创建jjVolcano火山图3（环形）...")

  pdf(file = file.path(output_dirs$scRNAtoolVis, "05_jjVolcano_polar.pdf"), width = 8, height = 8)
  print(
    jjVolcano(diffData = markers,
              legend.position = c(0.9, 0.85),
              polar = TRUE,
              size = 3,
              fontface = 'italic',
              tile.col = jj_colors,
              log2FC.cutoff = analysis_params$logFCfilter) +
      labs(x = "Cell type")
  )
  dev.off()

  message("jjVolcano plot 3 saved / jjVolcano火山图3已保存")

} else {
  message("Warning: markers data not available, skipping volcano plots")
  message("警告：markers数据不可用，跳过火山图绑制")
}

message("scRNAtoolVis visualization completed!")
message("scRNAtoolVis可视化完成！")

# ==============================================================================
# 18. Single Cell Gene Knockout Simulation Analysis (scTenifoldKnk)
# 单细胞基因敲除模拟分析 - 仅在疾病组/实验组进行分析
# ==============================================================================

message("\n================================================================================")
message("         Step 18: Single Cell Gene Knockout Simulation Analysis                 ")
message("         单细胞基因敲除模拟分析 - 仅在疾病组/实验组进行                          ")
message("================================================================================")

if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  message("警告：scTenifoldKnk 包未安装，跳过基因敲除模拟分析")
  message("安装方法：devtools::install_github('cailab-tamu/scTenifoldKnk')")
} else {

  library(scTenifoldKnk)
  library(ggrepel)
  set.seed(123)

  # =============================================
  # 18.1 Set knockout analysis parameters / 设置敲除分析参数
  # =============================================

  knockout_params <- list(
    knockout_genes = showGenes,
    n_hvg = 2000,
    qc_mtThreshold = 0.1,
    qc_minLSize = 1000,
    nc_nNet = 3,
    nc_nCells = 300,
    pval_threshold = 0.05
  )

  message(sprintf("敲除目标基因: %s", paste(knockout_params$knockout_genes, collapse = ", ")))
  message(sprintf("用于网络构建的高变基因数: %d", knockout_params$n_hvg))

  # =============================================
  # 18.2 Create output directory / 创建输出目录
  # =============================================

  knockout_output_dir <- file.path(workDir, output_dirs$gene_knockout)
  if (!dir.exists(knockout_output_dir)) {
    dir.create(knockout_output_dir, recursive = TRUE)
  }

  # =============================================
  # 18.3 提取疾病组细胞进行分析
  # =============================================

  message("提取疾病组细胞...")

  scObject_disease_ko <- subset(Peripheral_Blood_Mononuclear_Cells, subset = Type == "Disease")

  message(sprintf("疾病组细胞数: %d", ncol(scObject_disease_ko)))
  message(sprintf("疾病组细胞类型: %s", paste(unique(Idents(scObject_disease_ko)), collapse = ", ")))

  # =============================================
  # 18.4 对每个敲除基因进行分析
  # =============================================

  available_knockout_genes <- knockout_params$knockout_genes[
    knockout_params$knockout_genes %in% rownames(scObject_disease_ko)
  ]
  missing_knockout_genes <- setdiff(knockout_params$knockout_genes, available_knockout_genes)

  if (length(missing_knockout_genes) > 0) {
    message(sprintf("警告: 以下基因不在数据中: %s", paste(missing_knockout_genes, collapse = ", ")))
  }

  if (length(available_knockout_genes) > 0) {

    message(sprintf("开始对 %d 个基因进行敲除模拟分析（仅疾病组）", length(available_knockout_genes)))

    countMat <- GetAssayData(scObject_disease_ko, layer = "counts")

    scObject_disease_ko <- FindVariableFeatures(
      object = scObject_disease_ko, selection.method = "vst",
      nfeatures = knockout_params$n_hvg
    )
    hvgs <- VariableFeatures(scObject_disease_ko)

    for (target_gene in available_knockout_genes) {

      message(sprintf("\n########## 分析基因: %s ##########", target_gene))

      gene_output_dir <- file.path(knockout_output_dir, paste0("Gene_", target_gene))
      if (!dir.exists(gene_output_dir)) {
        dir.create(gene_output_dir, recursive = TRUE)
      }

      tryCatch({
        data <- as.data.frame(countMat[unique(c(target_gene, hvgs)), ])
        message(sprintf("  表达矩阵: %d 基因 x %d 细胞", nrow(data), ncol(data)))

        message("  运行scTenifoldKnk分析...")
        result <- scTenifoldKnk(
          countMatrix = data,
          gKO = target_gene,
          qc_mtThreshold = knockout_params$qc_mtThreshold,
          qc_minLSize = knockout_params$qc_minLSize,
          nc_nNet = knockout_params$nc_nNet,
          nc_nCells = knockout_params$nc_nCells
        )
        message("  scTenifoldKnk分析完成!")

        df <- result$diffRegulation
        df <- df[df$gene != target_gene, ]

        outTab <- df[df$p.adj < knockout_params$pval_threshold, ]
        write.csv(outTab, file = file.path(gene_output_dir, paste0(target_gene, "_sigDiff.csv")),
                  row.names = FALSE)
        message(sprintf("  发现 %d 个显著差异调控基因 (p.adj < 0.05)", nrow(outTab)))

        # =============================================
        # 散点图（Scatter Plot）
        # =============================================
        message("  生成散点图...")

        df$log_p.adj <- -log10(df$p.adj)
        df$significant <- ifelse(df$p.adj < knockout_params$pval_threshold, "Significant", "Not significant")
        label_genes <- subset(df, p.adj < knockout_params$pval_threshold)

        y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)

        scatter_colors <- colorRampPalette(brewer.pal(12, "Set3"))(2)
        p_scatter <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
          geom_point(alpha = 0.7, size = 1.5) +
          scale_color_manual(values = c("Significant" = scatter_colors[1], "Not significant" = "gray70")) +
          geom_hline(yintercept = -log10(knockout_params$pval_threshold), linetype = "dashed", color = scatter_colors[1]) +
          geom_text_repel(data = label_genes, aes(label = gene), size = 3, max.overlaps = 50,
                          color = "black", fontface = "italic") +
          labs(title = paste0(target_gene, " Knockout"),
               x = "Z-score",
               y = "-log10(p.adj)") +
          theme_classic(base_size = 14) +
          coord_cartesian(ylim = c(0, y_upper)) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12),
            legend.position = "none"
          )

        pdf(file.path(gene_output_dir, paste0(target_gene, "_KO_scatter.pdf")), width = 6, height = 5)
        print(p_scatter)
        dev.off()
        message("  已保存: ", target_gene, "_KO_scatter.pdf")

        # =============================================
        # 柱状图（Bar Plot）- Top 20基因
        # =============================================
        message("  生成柱状图...")

        top_genes <- head(df[order(-df$FC), ], 20)
        bar_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nrow(top_genes))

        p_bar <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC, fill = reorder(gene, FC))) +
          geom_bar(stat = 'identity', alpha = 0.9) +
          scale_fill_manual(values = bar_colors) +
          coord_flip() +
          labs(title = paste0("Top 20 Differentially Regulated Genes\n(", target_gene, " Knockout)"),
               x = "Gene",
               y = "FC") +
          theme_classic(base_size = 14) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.title = element_text(size = 13, face = "bold"),
            axis.text.y = element_text(size = 11, face = "italic"),
            axis.text.x = element_text(size = 11),
            legend.position = "none"
          )

        pdf(file.path(gene_output_dir, paste0(target_gene, "_KO_barplot.pdf")), width = 6, height = 5)
        print(p_bar)
        dev.off()
        message("  已保存: ", target_gene, "_KO_barplot.pdf")

        # 保存完整结果
        write.csv(df, file = file.path(gene_output_dir, paste0(target_gene, "_KO_allResults.csv")),
                  row.names = FALSE)

        message(sprintf("  基因 %s 敲除分析完成!", target_gene))

      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        message(sprintf("  基因 %s 敲除分析失败: %s", target_gene, e$message))
      })
    }

  } else {
    message("警告: 没有有效的敲除基因!")
  }

  message("\n================================================================================")
  message("         基因敲除模拟分析完成!                                                  ")
  message("================================================================================")
  message(sprintf("结果保存至: %s", knockout_output_dir))

  message("\n=== 敲除分析输出结构 ===")
  message("21.Gene_Knockout_Analysis/")
  message("  └── Gene_[GeneName]/")
  message("      ├── [Gene]_KO_scatter.pdf      (散点图)")
  message("      ├── [Gene]_KO_barplot.pdf      (柱状图)")
  message("      ├── [Gene]_sigDiff.csv         (显著差异基因)")
  message("      └── [Gene]_KO_allResults.csv   (完整结果)")
}

# ==============================================================================
# 17.5 Cell Type Marker Gene Correlation Bubble Plot
# 细胞类型标记基因相关性气泡图
# ==============================================================================

message("\n================================================================================")
message("         Cell Type Marker Gene Correlation Bubble Plot                          ")
message("         细胞类型标记基因相关性气泡图                                            ")
message("================================================================================")

# 创建输出目录
corr_bubble_dir <- file.path(workDir, "22.Correlation_BubblePlot")
if (!dir.exists(corr_bubble_dir)) dir.create(corr_bubble_dir, recursive = TRUE)

tryCatch({
  # 获取表达矩阵和细胞类型信息
  expr_mat <- GetAssayData(Peripheral_Blood_Mononuclear_Cells, layer = "data")
  cell_types <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
  unique_celltypes <- sort(unique(cell_types))

  # 使用FindAllMarkers的结果获取每个细胞类型的marker基因
  if (exists("sig_markers") && nrow(sig_markers) > 0) {
    marker_data <- sig_markers
  } else if (exists("markers") && nrow(markers) > 0) {
    marker_data <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0, ]
  } else {
    stop("No marker gene data available")
  }

  message(sprintf("共有 %d 个细胞类型, %d 个显著marker基因",
                  length(unique_celltypes), nrow(marker_data)))

  # 计算每个marker基因与其对应细胞类型的相关性
  # 对每个细胞类型，构建二元向量(是否属于该类型)，计算与基因表达的相关性
  top_n_genes <- 3  # 每个细胞类型选取相关性前3的基因

  corr_results <- data.frame()

  for (ct in unique_celltypes) {
    # 获取该细胞类型的marker基因
    ct_markers <- marker_data[marker_data$cluster == ct, ]
    if (nrow(ct_markers) == 0) next

    # 按avg_log2FC排序，先取候选基因（最多取前20个计算相关性）
    ct_markers <- ct_markers[order(-ct_markers$avg_log2FC), ]
    candidate_genes <- head(ct_markers$gene, 20)
    candidate_genes <- candidate_genes[candidate_genes %in% rownames(expr_mat)]
    if (length(candidate_genes) == 0) next

    # 构建该细胞类型的二元向量
    ct_binary <- as.numeric(cell_types == ct)

    # 计算每个候选基因与该细胞类型的相关性
    gene_corrs <- sapply(candidate_genes, function(g) {
      gene_expr <- as.numeric(expr_mat[g, ])
      cor(gene_expr, ct_binary, method = "spearman")
    })

    # 选取相关性最高的前N个基因
    gene_corrs <- sort(gene_corrs, decreasing = TRUE)
    selected_genes <- head(names(gene_corrs), top_n_genes)

    # 计算选中基因在所有细胞类型中的相关性和表达百分比
    for (g in selected_genes) {
      gene_expr <- as.numeric(expr_mat[g, ])
      for (ct2 in unique_celltypes) {
        ct2_binary <- as.numeric(cell_types == ct2)
        r_val <- cor(gene_expr, ct2_binary, method = "spearman")
        # 计算该基因在该细胞类型中的表达百分比
        ct2_cells <- which(cell_types == ct2)
        pct_expr <- sum(gene_expr[ct2_cells] > 0) / length(ct2_cells) * 100
        # 计算平均表达量
        avg_expr <- mean(gene_expr[ct2_cells])

        corr_results <- rbind(corr_results, data.frame(
          Gene = g,
          MarkerFor = ct,
          CellType = ct2,
          Correlation = r_val,
          Pct_Expressed = pct_expr,
          Avg_Expression = avg_expr,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  if (nrow(corr_results) == 0) stop("No correlation results computed")

  message(sprintf("计算完成: %d 个基因-细胞类型组合", nrow(corr_results)))

  # 保存相关性结果
  write.csv(corr_results, file.path(corr_bubble_dir, "01_Correlation_Results.csv"), row.names = FALSE)

  # 设置基因顺序：按所属细胞类型分组
  gene_order <- unique(corr_results[order(corr_results$MarkerFor, -corr_results$Correlation), "Gene"])
  corr_results$Gene <- factor(corr_results$Gene, levels = rev(gene_order))
  corr_results$CellType <- factor(corr_results$CellType, levels = unique_celltypes)

  # 为每个基因添加所属细胞类型标签（用于分面或颜色标注）
  gene_to_ct <- corr_results[!duplicated(corr_results$Gene), c("Gene", "MarkerFor")]

  # 绘制相关性气泡图
  n_genes <- length(unique(corr_results$Gene))
  n_cts <- length(unique_celltypes)
  plot_height <- max(8, n_genes * 0.35 + 2)
  plot_width <- max(8, n_cts * 0.8 + 3)

  # 生成细胞类型对应的颜色（用于左侧基因分组标注）
  ct_colors <- colorRampPalette(brewer.pal(min(12, length(unique_celltypes)), "Set3"))(length(unique_celltypes))
  names(ct_colors) <- unique_celltypes

  p_bubble <- ggplot(corr_results, aes(x = CellType, y = Gene)) +
    geom_point(aes(size = Pct_Expressed, color = Correlation)) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0, limits = c(-0.5, 0.8),
                          oob = scales::squish,
                          name = "Correlation") +
    scale_size_continuous(range = c(0.5, 6), name = "% Expressed",
                          breaks = c(10, 25, 50, 75, 100)) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 9, face = "italic"),
      axis.title = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = NULL)

  # 添加左侧细胞类型分组色条标注
  # 创建基因-细胞类型映射用于左侧标注
  gene_ct_map <- gene_to_ct
  gene_ct_map$Gene <- factor(gene_ct_map$Gene, levels = rev(gene_order))
  gene_ct_map <- gene_ct_map[order(gene_ct_map$Gene), ]

  p_strip <- ggplot(gene_ct_map, aes(x = 1, y = Gene, fill = MarkerFor)) +
    geom_tile(width = 0.8) +
    scale_fill_manual(values = ct_colors, name = "Marker For") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )

  # 组合图形
  p_combined <- p_strip + p_bubble +
    plot_layout(widths = c(0.03, 1)) +
    plot_annotation(title = "Cell Type Marker Gene Correlation Bubble Plot",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))

  pdf(file.path(corr_bubble_dir, "02_CellType_Marker_Correlation_BubblePlot.pdf"),
      width = plot_width, height = plot_height)
  print(p_combined)
  dev.off()

  message("相关性气泡图已保存至: 22.Correlation_BubblePlot/")
  message("Correlation bubble plot saved!")

  # ---- 分组绘制：对照组和疾病组分别绘制 ----
  message("绘制分组相关性气泡图（对照组 vs 疾病组）...")

  group_info <- Peripheral_Blood_Mononuclear_Cells$Type

  for (grp in c("control", "Disease")) {
    grp_cells <- which(group_info == grp)
    if (length(grp_cells) < 10) {
      message(sprintf("  %s 组细胞数不足，跳过", grp))
      next
    }

    grp_expr <- expr_mat[, grp_cells]
    grp_celltypes <- cell_types[grp_cells]
    grp_unique_cts <- sort(unique(grp_celltypes))

    grp_corr <- data.frame()

    for (g in gene_order) {
      if (!g %in% rownames(grp_expr)) next
      gene_expr <- as.numeric(grp_expr[g, ])
      marker_ct <- gene_to_ct$MarkerFor[gene_to_ct$Gene == g]

      for (ct2 in grp_unique_cts) {
        ct2_binary <- as.numeric(grp_celltypes == ct2)
        r_val <- cor(gene_expr, ct2_binary, method = "spearman")
        ct2_idx <- which(grp_celltypes == ct2)
        pct_expr <- sum(gene_expr[ct2_idx] > 0) / length(ct2_idx) * 100

        grp_corr <- rbind(grp_corr, data.frame(
          Gene = g, MarkerFor = marker_ct, CellType = ct2,
          Correlation = r_val, Pct_Expressed = pct_expr,
          stringsAsFactors = FALSE
        ))
      }
    }

    if (nrow(grp_corr) == 0) next

    grp_corr$Gene <- factor(grp_corr$Gene, levels = rev(gene_order))
    grp_corr$CellType <- factor(grp_corr$CellType, levels = grp_unique_cts)

    p_grp <- ggplot(grp_corr, aes(x = CellType, y = Gene)) +
      geom_point(aes(size = Pct_Expressed, color = Correlation)) +
      scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                            midpoint = 0, limits = c(-0.5, 0.8),
                            oob = scales::squish, name = "Correlation") +
      scale_size_continuous(range = c(0.5, 6), name = "% Expressed",
                            breaks = c(10, 25, 50, 75, 100)) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 9, face = "italic"),
        axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      ) +
      labs(title = paste0("Marker Gene Correlation - ", grp, " Group"))

    grp_height <- max(8, length(gene_order) * 0.35 + 2)
    grp_width <- max(8, length(grp_unique_cts) * 0.8 + 3)

    pdf(file.path(corr_bubble_dir, paste0("03_Correlation_BubblePlot_", grp, ".pdf")),
        width = grp_width, height = grp_height)
    print(p_grp)
    dev.off()

    message(sprintf("  %s 组相关性气泡图已保存", grp))
  }

  message("分组相关性气泡图绑制完成!")

}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  message(sprintf("相关性气泡图绑制失败: %s", e$message))
})

# ==============================================================================
# 17.6 ShowGenes Module Score Violin Plot (Disease Group Only)
# 研究基因集整体评分在疾病组各细胞类型中的小提琴图
# ==============================================================================

message("\n================================================================================")
message("         ShowGenes Module Score Violin Plot (Disease Group Only)                ")
message("         研究基因集整体评分 - 疾病组各细胞类型小提琴图                            ")
message("================================================================================")

violin_output_dir <- file.path(workDir, "23.ShowGenes_VlnPlot")
if (!dir.exists(violin_output_dir)) dir.create(violin_output_dir, recursive = TRUE)

tryCatch({
  # 确保Idents为细胞类型
  Idents(Peripheral_Blood_Mononuclear_Cells) <- Peripheral_Blood_Mononuclear_Cells$cell_type

  # 筛选在数据中存在的showGenes
  valid_genes <- showGenes[showGenes %in% rownames(Peripheral_Blood_Mononuclear_Cells)]
  if (length(valid_genes) == 0) stop("showGenes中没有基因存在于表达矩阵中")
  message(sprintf("有效基因: %d / %d (%s)", length(valid_genes), length(showGenes),
                  paste(valid_genes, collapse = ", ")))

  # 计算基因集整体Module Score
  Peripheral_Blood_Mononuclear_Cells <- AddModuleScore(
    Peripheral_Blood_Mononuclear_Cells,
    features = list(valid_genes),
    name = "ShowGenes_Score"
  )
  message("Module Score计算完成")

  # 只提取疾病组细胞
  disease_obj <- subset(Peripheral_Blood_Mononuclear_Cells, subset = Type == "Disease")
  Idents(disease_obj) <- disease_obj$cell_type
  n_disease <- ncol(disease_obj)
  message(sprintf("疾病组细胞数: %d", n_disease))

  unique_cts <- sort(unique(as.character(Idents(disease_obj))))
  n_cts <- length(unique_cts)

  # 构建绑图数据
  plot_df <- data.frame(
    CellType = as.character(Idents(disease_obj)),
    ModuleScore = disease_obj$ShowGenes_Score1,
    stringsAsFactors = FALSE
  )
  plot_df$CellType <- factor(plot_df$CellType, levels = unique_cts)

  # 生成每个细胞类型不同颜色
  if (n_cts <= 9) {
    ct_cols <- brewer.pal(max(3, n_cts), "Set1")[1:n_cts]
  } else {
    ct_cols <- colorRampPalette(brewer.pal(9, "Set1"))(n_cts)
  }
  names(ct_cols) <- unique_cts

  # 绘制小提琴图
  p_module_vln <- ggplot(plot_df, aes(x = CellType, y = ModuleScore, fill = CellType)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.85, linewidth = 0.4) +
    geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.5,
                 outlier.alpha = 0.4, linewidth = 0.4) +
    scale_fill_manual(values = ct_cols) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40")
    ) +
    labs(y = "Module Score",
         title = "ShowGenes Module Score (Disease Group)",
         subtitle = paste0("Genes: ", paste(valid_genes, collapse = ", ")))

  vln_width <- max(7, n_cts * 0.9 + 2)

  pdf(file.path(violin_output_dir, "01_ShowGenes_ModuleScore_VlnPlot_Disease.pdf"),
      width = vln_width, height = 6)
  print(p_module_vln)
  dev.off()

  # 保存数据
  write.csv(plot_df, file.path(violin_output_dir, "01_ModuleScore_Data_Disease.csv"),
            row.names = FALSE)

  message("疾病组基因集整体评分小提琴图已保存!")

  # ---- 单个基因小提琴图（疾病组）----
  message("绑制单个基因小提琴图（疾病组）...")

  indiv_dir <- file.path(violin_output_dir, "Individual_Genes_Disease")
  if (!dir.exists(indiv_dir)) dir.create(indiv_dir, recursive = TRUE)

  expr_disease <- GetAssayData(disease_obj, layer = "data")

  for (g in valid_genes) {
    tryCatch({
      gene_expr <- as.numeric(expr_disease[g, ])
      g_df <- data.frame(
        CellType = factor(as.character(Idents(disease_obj)), levels = unique_cts),
        Expression = gene_expr,
        stringsAsFactors = FALSE
      )

      p_g <- ggplot(g_df, aes(x = CellType, y = Expression, fill = CellType)) +
        geom_violin(trim = FALSE, scale = "width", alpha = 0.85, linewidth = 0.4) +
        geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.5,
                     outlier.alpha = 0.4, linewidth = 0.4) +
        scale_fill_manual(values = ct_cols) +
        theme_classic(base_size = 13) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14)
        ) +
        labs(y = "Expression Level", title = g)

      pdf(file.path(indiv_dir, paste0(g, "_VlnPlot_Disease.pdf")),
          width = max(7, n_cts * 0.9 + 2), height = 5.5)
      print(p_g)
      dev.off()

    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message(sprintf("  基因 %s 小提琴图失败: %s", g, e$message))
    })
  }

  message(sprintf("单个基因小提琴图已保存至: %s", indiv_dir))

}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  message(sprintf("ShowGenes Module Score小提琴图绑制失败: %s", e$message))
})

# ==============================================================================
# 17.7 Ro/e (Observed/Expected Ratio) Heatmap for Cell Types
# 细胞类型Ro/e观察/期望比值热图
# ==============================================================================

message("\n================================================================================")
message("         Ro/e (Observed/Expected Ratio) Heatmap                                 ")
message("         细胞类型Ro/e观察/期望比值热图                                           ")
message("================================================================================")

roe_output_dir <- file.path(workDir, "24.Roe_Heatmap")
if (!dir.exists(roe_output_dir)) dir.create(roe_output_dir, recursive = TRUE)

tryCatch({
  # 确保Idents为细胞类型
  Idents(Peripheral_Blood_Mononuclear_Cells) <- Peripheral_Blood_Mononuclear_Cells$cell_type

  cell_types <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
  groups <- Peripheral_Blood_Mononuclear_Cells$Type

  # 构建观察频数表: 行=细胞类型, 列=分组
  obs_table <- table(cell_types, groups)
  message("观察频数表:")
  print(obs_table)

  # 计算期望频数 (基于卡方检验的期望值)
  # E_ij = (行合计_i * 列合计_j) / 总数
  row_sums <- rowSums(obs_table)
  col_sums <- colSums(obs_table)
  total <- sum(obs_table)

  exp_table <- outer(row_sums, col_sums) / total

  # 计算 Ro/e
  roe_table <- obs_table / exp_table

  message("\nRo/e比值表:")
  print(round(roe_table, 3))

  # 判断富集/缺失标记: Ro/e > 1 为 "+", Ro/e ≈ 1 为 "±", Ro/e < 1 为 "-"
  label_table <- matrix("", nrow = nrow(roe_table), ncol = ncol(roe_table))
  for (i in 1:nrow(roe_table)) {
    for (j in 1:ncol(roe_table)) {
      val <- roe_table[i, j]
      if (val > 1.1) {
        label_table[i, j] <- "+"
      } else if (val < 0.9) {
        label_table[i, j] <- "-"
      } else {
        label_table[i, j] <- "±"
      }
    }
  }

  # 保存Ro/e结果
  roe_df <- as.data.frame.matrix(roe_table)
  roe_df$CellType <- rownames(roe_df)
  roe_df <- roe_df[, c("CellType", colnames(obs_table))]
  write.csv(roe_df, file.path(roe_output_dir, "01_Roe_Values.csv"), row.names = FALSE)

  # 使用ggplot2绘制Ro/e热图（避免ComplexHeatmap的PDF兼容性问题）
  roe_mat <- as.matrix(roe_table)

  # 构建长格式数据用于ggplot
  roe_long <- data.frame(
    CellType = rep(rownames(roe_mat), ncol(roe_mat)),
    Group = rep(colnames(roe_mat), each = nrow(roe_mat)),
    Roe = as.vector(roe_mat),
    stringsAsFactors = FALSE
  )

  # 添加富集标记
  roe_long$Label <- ifelse(roe_long$Roe > 1.1, "+",
                    ifelse(roe_long$Roe < 0.9, "-", "\u00b1"))
  roe_long$DisplayText <- paste0(round(roe_long$Roe, 2), "\n", roe_long$Label)

  # 设置因子顺序
  roe_long$CellType <- factor(roe_long$CellType, levels = rev(sort(unique(roe_long$CellType))))
  roe_long$Group <- factor(roe_long$Group, levels = c("control", "Disease"))

  n_ct <- length(unique(roe_long$CellType))
  ht_height <- max(5, n_ct * 0.6 + 2)

  p_roe <- ggplot(roe_long, aes(x = Group, y = CellType, fill = Roe)) +
    geom_tile(color = "grey80", linewidth = 0.5) +
    geom_text(aes(label = DisplayText,
                  color = ifelse(Roe > 2, "white", "black")),
              size = 3.5, lineheight = 0.85) +
    scale_color_identity() +
    scale_fill_gradient2(low = "#F7FBFF", mid = "#9ECAE1", high = "#08519C",
                         midpoint = 1.5, limits = c(0, 3),
                         oob = scales::squish, name = "Ro/e") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    labs(title = "Ro/e: Observed/Expected Cell Ratio (Control vs Disease)")

  pdf(file.path(roe_output_dir, "02_Roe_Heatmap.pdf"),
      width = 6, height = ht_height)
  print(p_roe)
  dev.off()

  message("Ro/e热图已保存至: 24.Roe_Heatmap/")
  message("Ro/e heatmap saved!")

}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  message(sprintf("Ro/e热图绑制失败: %s", e$message))
})

# ==============================================================================
# 17. Statistical Report Generation
# ==============================================================================

message("Generating statistical reports")

# Create comprehensive statistics table
stats_summary <- data.frame(
  Parameter = c(
    "Total_Cells_Raw", "Total_Cells_Filtered",
    "Total_Genes_Raw", "Total_Genes_Filtered",
    "Variable_Genes", "Cell_Types_Identified",
    "control_Cells", "Disease_Cells",
    "PC_Used", "Clustering_Resolution"
  ),
  Value = c(
    cells_before, cells_after,
    genes_before, genes_after,
    length(VariableFeatures(Peripheral_Blood_Mononuclear_Cells)),
    length(unique(Idents(Peripheral_Blood_Mononuclear_Cells))),
    sum(Peripheral_Blood_Mononuclear_Cells$Type == "control"),
    sum(Peripheral_Blood_Mononuclear_Cells$Type == "Disease"),
    analysis_params$pcSelect,
    analysis_params$cluster_resolution
  )
)

write.csv(stats_summary,
          file.path(output_dirs$statistics_reports, "01_Analysis_Summary.csv"),
          row.names = FALSE)

# ==============================================================================
# 15. Save Final Results
# ==============================================================================

message("Saving final analysis results")

# Save Seurat object
cellAnn <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
names(cellAnn) <- colnames(Peripheral_Blood_Mononuclear_Cells)

save(Peripheral_Blood_Mononuclear_Cells, cellAnn,
     file = file.path(output_dirs$final_results, "01_SingleCell_Analysis_Complete.RData"))

# Save filtered expression matrix
filtered_counts_matrix <- GetAssayData(Peripheral_Blood_Mononuclear_Cells, layer = "counts")
saveRDS(filtered_counts_matrix,
        file.path(output_dirs$final_results, "02_Filtered_Expression_Matrix.rds"))

# Save cell annotation mapping table
cell_mapping <- data.frame(
  Cell_ID = colnames(Peripheral_Blood_Mononuclear_Cells),
  Cluster_ID = as.character(Peripheral_Blood_Mononuclear_Cells$seurat_clusters),
  Cell_Type = as.character(Idents(Peripheral_Blood_Mononuclear_Cells)),
  Group = Peripheral_Blood_Mononuclear_Cells$Type,
  nFeature_RNA = Peripheral_Blood_Mononuclear_Cells$nFeature_RNA,
  nCount_RNA = Peripheral_Blood_Mononuclear_Cells$nCount_RNA,
  percent_mt = Peripheral_Blood_Mononuclear_Cells$percent.mt,
  stringsAsFactors = FALSE
)

write.csv(cell_mapping,
          file.path(output_dirs$final_results, "03_Cell_ID_to_CellType_Mapping.csv"),
          row.names = FALSE)

# ==============================================================================
# 16. Analysis Completion Report
# ==============================================================================

message("\n================================================================================")
message("                  Single Cell RNA-seq Comprehensive Analysis Complete!          ")
message("                      单细胞RNA测序综合分析完成！                                ")

