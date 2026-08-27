
# ==================================================================================
# 17_wgcna_analysis.R
# WGCNA (Weighted Gene Co-expression Network Analysis) on RCC transcriptomic
# data, identifying modules associated with consensus clustering subtypes and
# ferro-aging-related candidate core genes.
# ==================================================================================
# ------------------------ Package loading ---------------------------
cat("==== Loading R packages ====\n")
library(WGCNA)
library(limma)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(scales)
library(viridis)
library(corrplot)
library(cowplot)

# 显式指定使用stats包的dist函数，避免冲突
dist <- stats::dist

# ==================================================================================
# ------------------------ 参数设置 ---------------------------
# Set working directory and input file paths (relative paths only; no absolute paths used)
workDir <- getwd()
expFilePath <- "merged_after_batch_removal.csv"
clusterFilePath <- "geneCluster.txt"  # subtype file path

# 分析参数
min_module_size <- 50                    # 最小模块大小 (建议30-100)
merge_threshold <- 0.25                  # 模块合并阈值 (建议0.15-0.30)
soft_threshold_range <- 1:30             # 软阈值范围
cutHeight_sample <- 15000                # 样本聚类切割高度

# 输出参数
figure_width <- 12                       # 图片宽度 (增加宽度)
figure_height <- 8                       # 图片高度
figure_dpi <- 300                        # 图片分辨率
pdf_pointsize <- 12                      # PDF字体大小

# 统计显著性阈值
p_threshold_high <- 0.001                # 高显著性 (***)
p_threshold_medium <- 0.01               # 中等显著性 (**)
p_threshold_low <- 0.05                  # 低显著性 (*)

# 颜色设置
color_palette <- "RdBu"                  # 热图配色方案
use_custom_colors <- TRUE                # 是否使用自定义颜色

# ==================================================================================
# ------------------------ 环境初始化 ---------------------------
cat("==== WGCNA高级版分析开始 ====\n")
start_time <- Sys.time()

# 保存参数到变量（避免清理环境时丢失）
temp_params <- list(
  workDir = workDir,
  expFilePath = expFilePath,
  clusterFilePath = clusterFilePath,
  min_module_size = min_module_size,
  merge_threshold = merge_threshold,
  soft_threshold_range = soft_threshold_range,
  cutHeight_sample = cutHeight_sample,
  figure_width = figure_width,
  figure_height = figure_height,
  figure_dpi = figure_dpi,
  pdf_pointsize = pdf_pointsize,
  p_threshold_high = p_threshold_high,
  p_threshold_medium = p_threshold_medium,
  p_threshold_low = p_threshold_low,
  color_palette = color_palette,
  use_custom_colors = use_custom_colors,
  start_time = start_time  # 保存开始时间
)

# 清理环境
rm(list = setdiff(ls(), "temp_params"))
gc()

# 变量
list2env(temp_params, envir = .GlobalEnv)
rm(temp_params)

# 加载包
suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(scales)
  library(viridis)
  library(corrplot)
  library(cowplot)
})

# 显式指定使用stats包的dist函数，避免冲突
dist <- stats::dist

# 创建输出目录结构
output_dirs <- c(
  "01_Quality_Control",          # 质量控制图表
  "02_Network_Analysis",         # 网络分析图表
  "03_Module_Detection",         # 模块检测图表
  "04_Module_Trait_Analysis",    # 模块-性状关联分析
  "05_Gene_Analysis",           # 基因水平分析
  "06_Data_Output",             # 数据输出文件
  "07_Module_Visualization"     # 模块可视化
)

# 确保在工作目录中创建目录
current_wd <- getwd()
setwd(workDir)  # 切换到工作目录

for(dir in output_dirs) {
  if(!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("创建目录:", dir, "\n")
  }
}

cat("当前工作目录:", getwd(), "\n")
cat("输出目录创建完成\n")

# 主题设置
theme_sci <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 1),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold")
    )
}

# 设置随机种子和进度跟踪
set.seed(12345)

# 进度跟踪
total_steps <- 12
current_step <- 0
update_progress <- function(message) {
  current_step <<- current_step + 1
  cat(sprintf("[%d/%d] %s\n", current_step, total_steps, message))
}

# ------------------------ 数据读取与预处理 ---------------------------
update_progress("数据读取与预处理")

# 检查文件是否存在
if (!file.exists(expFilePath)) {
  stop("表达矩阵文件不存在！请检查路径: ", file.path(getwd(), expFilePath))
}
if (!file.exists(clusterFilePath)) {
  stop("分型文件不存在！请检查路径: ", file.path(getwd(), clusterFilePath))
}

# 读取分型文件
cat("读取分型文件...\n")
clusterData <- read.table(clusterFilePath, sep = "\t", header = TRUE,
                          check.names = FALSE, stringsAsFactors = FALSE)
colnames(clusterData) <- c("ID", "ClusterType")
cat("分型文件包含", nrow(clusterData), "个样本\n")

# 检测分型数量
unique_clusters <- unique(clusterData$ClusterType)
num_clusters <- length(unique_clusters)
cat("检测到", num_clusters, "个分型:", paste(unique_clusters, collapse = ", "), "\n")

# 验证分型数量在2-4之间
if (num_clusters < 2 || num_clusters > 4) {
  stop("分型数量必须在2-4之间，当前检测到", num_clusters, "个分型")
}

# 统计每个分型的样本数
for (cluster in unique_clusters) {
  count <- sum(clusterData$ClusterType == cluster)
  cat("分型", cluster, "样本数:", count, "\n")
}

# 自动检测分隔符并读取表达数据
tmp_head <- readLines(expFilePath, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")
rawData <- read.table(expFilePath, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(rawData) <- rawData[, 1]
exprData <- rawData[, -1]

# 数据预处理流程
dimNames <- list(rownames(exprData), colnames(exprData))
numericData <- matrix(as.numeric(as.matrix(exprData)), nrow = nrow(exprData), dimnames = dimNames)

# Log2转换
if (max(numericData, na.rm = TRUE) > 1000) {
  numericData <- log2(numericData + 1)
  cat("已进行Log2转换\n")
} else {
  cat("数据似乎已经过Log2转换，跳过\n")
}

# limma归一化
normalizedData <- normalizeBetweenArrays(numericData)

# 过滤低方差基因
gene_sd <- apply(normalizedData, 1, sd, na.rm = TRUE)
filteredData <- normalizedData[gene_sd > 0.25, ]
cat("过滤后保留基因数:", nrow(filteredData), "\n")

# 匹配表达数据和分型数据
sample_names <- colnames(filteredData)
matched_samples <- intersect(sample_names, clusterData$ID)
cat("表达数据和分型数据共有", length(matched_samples), "个样本\n")

if (length(matched_samples) == 0) {
  stop("表达数据和分型数据没有共同样本！请检查样本ID是否一致")
}

# 只保留有分型信息的样本
filteredData <- filteredData[, matched_samples]
clusterData <- clusterData[clusterData$ID %in% matched_samples, ]

# 确保样本顺序一致
clusterData <- clusterData[match(colnames(filteredData), clusterData$ID), ]

# 创建分组信息
group_info <- clusterData$ClusterType
names(group_info) <- clusterData$ID

cat("最终分析样本数:", ncol(filteredData), "\n")

# 数据转置为WGCNA格式
exprMatrix <- t(filteredData)

# ------------------------ 样本质量控制 ---------------------------
update_progress("样本质量控制")

# 样本和基因质量检查
enableWGCNAThreads()
sampleCheck <- goodSamplesGenes(exprMatrix, verbose = 3)
if (!sampleCheck$allOK) {
  exprMatrix <- exprMatrix[sampleCheck$goodSamples, sampleCheck$goodGenes]
  cat("移除低质量样本/基因后，维度:", dim(exprMatrix), "\n")
}

# 样本聚类分析
sampleDist <- dist(exprMatrix)
sampleDendro <- hclust(sampleDist, method = "average")

# 高质量样本聚类图（不显示样本名）
pdf(file.path("01_Quality_Control", "01_Sample_Clustering_Analysis.pdf"),
    width = figure_width, height = figure_height, pointsize = pdf_pointsize)
par(cex = 0.8, mar = c(4, 6, 4, 2))
plot(sampleDendro,
     main = "Sample Clustering for Outlier Detection",
     sub = "", xlab = "Samples", ylab = "Distance",
     cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.5,
     labels = FALSE)

# 为不同分型分配颜色
cluster_colors_map <- c("#3498DB", "#E74C3C", "#2ECC71", "#F39C12")
names(cluster_colors_map) <- unique_clusters[1:min(4, length(unique_clusters))]

# 获取每个样本的颜色
sample_colors <- cluster_colors_map[group_info[rownames(exprMatrix)]]

# 在底部添加颜色条
rect(xleft = 1:nrow(exprMatrix) - 0.4,
     xright = 1:nrow(exprMatrix) + 0.4,
     ybottom = -max(sampleDendro$height)*0.02,
     ytop = -max(sampleDendro$height)*0.01,
     col = sample_colors, border = sample_colors, xpd = TRUE)

# 添加图例
legend("topright", legend = unique_clusters,
       col = cluster_colors_map[unique_clusters], pch = 15, cex = 1.2)
dev.off()

# 异常样本检测
cutHeight <- cutHeight_sample
clusterCut <- cutreeStatic(sampleDendro, cutHeight = cutHeight, minSize = 3)
if (length(unique(clusterCut)) > 1) {
  keepSamples <- (clusterCut == 1)
  exprMatrix <- exprMatrix[keepSamples, ]
  cat("移除异常样本后，样本数:", nrow(exprMatrix), "\n")
}

# ------------------------ 临床性状数据准备 ---------------------------
update_progress("临床性状数据准备")

# 创建临床数据 - 为每个分型创建一列
clinicalData <- data.frame(row.names = rownames(exprMatrix))

for (cluster in unique_clusters) {
  clinicalData[[cluster]] <- ifelse(group_info[rownames(exprMatrix)] == cluster, 1, 0)
}

cat("临床数据维度:", dim(clinicalData), "\n")
cat("临床数据列名:", colnames(clinicalData), "\n")

# 样本-性状热图
sampleColors <- numbers2colors(clinicalData, signed = FALSE)
pdf(file.path("01_Quality_Control", "02_Sample_Trait_Heatmap.pdf"),
    width = figure_width + 2, height = figure_height)
plotDendroAndColors(sampleDendro, sampleColors,
                    groupLabels = names(clinicalData),
                    main = "Sample Dendrogram and Clinical Traits")
dev.off()

# ------------------------ 软阈值选择 ---------------------------
update_progress("软阈值选择与网络拓扑分析")

# 软阈值计算 (与原版一致: 默认unsigned网络)
powerVector <- soft_threshold_range
sftResult <- pickSoftThreshold(exprMatrix, powerVector = powerVector, verbose = 5)

if (is.null(sftResult$powerEstimate)) {
  optimalPower <- 6  # 默认值
  warning("无法自动确定最佳软阈值，使用默认值6")
} else {
  optimalPower <- sftResult$powerEstimate
}
cat("选择的软阈值:", optimalPower, "\n")

# 的软阈值选择图
sft_data <- data.frame(
  Power = sftResult$fitIndices[, 1],
  SFT_R2 = -sign(sftResult$fitIndices[, 3]) * sftResult$fitIndices[, 2],
  Mean_Connectivity = sftResult$fitIndices[, 5]
)

# 尺度自由拓扑图
p1 <- ggplot(sft_data, aes(x = Power, y = SFT_R2)) +
  geom_point(size = 3, color = "#2C3E50") +
  geom_text(aes(label = Power), hjust = -0.3, vjust = -0.3, size = 3.5) +
  geom_hline(yintercept = 0.9, color = "#E74C3C", linetype = "dashed", size = 1) +
  labs(x = "Soft Threshold (power)",
       y = "Scale Free Topology Model Fit (signed R²)",
       title = "Scale Independence") +
  theme_sci() +
  ylim(c(0, 1))

# 平均连通性图
p2 <- ggplot(sft_data, aes(x = Power, y = Mean_Connectivity)) +
  geom_point(size = 3, color = "#8E44AD") +
  geom_text(aes(label = Power), hjust = -0.3, vjust = -0.3, size = 3.5) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  theme_sci()

# 合并图表
pdf(file.path("02_Network_Analysis", "01_Soft_Threshold_Analysis.pdf"),
    width = figure_width + 2, height = figure_height - 2)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# ------------------------ 网络构建与模块检测 ---------------------------
update_progress("网络构建与TOM计算")

# 构建邻接矩阵和TOM 
adjacencyMatrix <- adjacency(exprMatrix, power = optimalPower)
TOMMatrix <- TOMsimilarity(adjacencyMatrix)
dissTOM <- 1 - TOMMatrix

# 基因聚类
geneDendro <- hclust(as.dist(dissTOM), method = "average")

# 基因聚类树图
pdf(file.path("02_Network_Analysis", "02_Gene_Clustering_Dendrogram.pdf"),
    width = figure_width + 5, height = figure_height)
par(cex = 0.8, mar = c(2, 6, 4, 2))
plot(geneDendro, xlab = "", sub = "",
     main = "Gene Clustering Based on TOM Dissimilarity",
     labels = FALSE, hang = 0.04,
     cex.main = 1.5, cex.lab = 1.3)
dev.off()

# ------------------------ 动态模块检测 ---------------------------
update_progress("动态模块检测与合并")

# 动态切树
minModuleSize <- min_module_size
dynamicMods <- cutreeDynamic(dendro = geneDendro, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
moduleColors <- labels2colors(dynamicMods)
cat("检测到", length(unique(moduleColors)), "个模块\n")

# 计算模块特征基因
MEs0 <- moduleEigengenes(exprMatrix, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# 模块合并
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
mergeThreshold <- merge_threshold

# 模块合并分析图
pdf(file.path("03_Module_Detection", "01_Module_Merging_Analysis.pdf"),
    width = figure_width, height = figure_height - 2)
par(cex = 0.8, mar = c(4, 6, 4, 2))
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "", cex.main = 1.5)
abline(h = mergeThreshold, col = "#E74C3C", lwd = 2, lty = 2)
text(x = length(METree$order)/2, y = mergeThreshold + 0.05,
     labels = paste("Merge threshold =", mergeThreshold),
     col = "#E74C3C", cex = 1.2)
dev.off()

# 执行模块合并
merge <- mergeCloseModules(exprMatrix, moduleColors, cutHeight = mergeThreshold, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# 合并前后对比图
pdf(file.path("03_Module_Detection", "02_Module_Colors_Comparison.pdf"),
    width = figure_width + 5, height = figure_height)
plotDendroAndColors(geneDendro, cbind(moduleColors, mergedColors),
                    c("Original", "Merged"), dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Module Assignment: Before and After Merging")
dev.off()

# 更新模块信息
moduleColors <- mergedColors
MEs <- mergedMEs
cat("合并后模块数:", length(unique(moduleColors)), "\n")

# ------------------------ 模块-性状关联分析 (经典WGCNA风格) ---------------------------
update_progress("模块-性状关联分析")

# 计算相关性
moduleTraitCor <- cor(MEs, clinicalData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(exprMatrix))

# 方法1: 经典WGCNA labeledHeatmap风格
pdf(file.path("04_Module_Trait_Analysis", "01_Module_Trait_Classic_WGCNA.pdf"),
    width = 8, height = 10)

# 设置图像参数
par(mar = c(6, 8.5, 3, 3))

# 准备文本矩阵（相关性系数和p值）
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# 创建经典WGCNA热图
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinicalData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               cex.lab.x = 1.2,
               cex.lab.y = 1.0,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# 方法2: 增强版经典风格（带模块颜色条）
pdf(file.path("04_Module_Trait_Analysis", "02_Module_Trait_Enhanced_Classic.pdf"),
    width = 10, height = 12)

# 准备数据
cor_matrix <- as.matrix(moduleTraitCor)
p_matrix <- as.matrix(moduleTraitPvalue)

# 创建文本标签（与参考图片一致的格式）
text_labels <- matrix(paste(sprintf("%.2f", cor_matrix),
                           "\n(",
                           ifelse(p_matrix < 1e-100, "< 1e-100",
                                 sprintf("%.0e", p_matrix)),
                           ")", sep = ""),
                     nrow = nrow(cor_matrix))

# 提取模块颜色
module_names <- rownames(cor_matrix)
module_colors <- gsub("ME", "", module_names)

# 创建模块颜色注释
color_annotation <- data.frame(
  Module = factor(module_colors, levels = module_colors)
)
rownames(color_annotation) <- module_names

# 创建颜色映射（使用实际的模块颜色）
module_color_map <- structure(module_colors, names = module_colors)

# 使用pheatmap创建增强热图
pheatmap(cor_matrix,
         # 颜色设置（经典蓝白红）
         color = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                                   "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                                   "#D6604D", "#B2182B", "#67001F"))(100),
         breaks = seq(-1, 1, length.out = 101),

         # 聚类设置
         cluster_rows = FALSE,
         cluster_cols = FALSE,

         # 显示设置
         show_rownames = TRUE,
         show_colnames = TRUE,

         # 字体设置
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 14,

         # 标题
         main = "Module-trait relationships",

         # 边框
         border_color = "white",

         # 单元格尺寸
         cellwidth = 80,
         cellheight = 35,

         # 显示数值
         display_numbers = text_labels,
         number_color = "black",
         fontsize_number = 10,

         # 行注释（模块颜色条）
         annotation_row = color_annotation,
         annotation_colors = list(Module = module_color_map),
         annotation_names_row = FALSE,

         # 图例设置
         legend = TRUE,
         legend_breaks = c(-1, -0.5, 0, 0.5, 1),
         legend_labels = c("-1", "-0.5", "0", "0.5", "1")
)
dev.off()

# 方法3: ggplot2版本的经典风格
library(ggplot2)
library(reshape2)

# 数据准备
cor_melted <- melt(moduleTraitCor)
p_melted <- melt(moduleTraitPvalue)
combined_data <- merge(cor_melted, p_melted, by = c("Var1", "Var2"))
colnames(combined_data) <- c("Module", "Trait", "Correlation", "Pvalue")

# 添加模块颜色
combined_data$ModuleColor <- gsub("ME", "", combined_data$Module)

# 创建文本标签（与参考图一致）
combined_data$Label <- paste(sprintf("%.2f", combined_data$Correlation),
                            "\n(",
                            ifelse(combined_data$Pvalue < 1e-100, "< 1e-100",
                                  sprintf("%.0e", combined_data$Pvalue)),
                            ")", sep = "")

# 确保模块顺序（与WGCNA输出一致）
module_order <- names(MEs)
combined_data$Module <- factor(combined_data$Module, levels = rev(module_order))

# 创建ggplot2版本的经典热图
p_classic <- ggplot(combined_data, aes(x = Trait, y = Module, fill = Correlation)) +
  # 基础热图瓦片
  geom_tile(color = "white", size = 0.8) +

  # 经典蓝白红配色方案
  scale_fill_gradient2(
    low = "#053061", mid = "white", high = "#67001F",
    midpoint = 0, limits = c(-1, 1),
    name = "",
    breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1")
  ) +

  # 添加数值标签
  geom_text(aes(label = Label), size = 3.2, color = "black", fontface = "bold") +

  # 主题设置
  theme_minimal() +
  theme(
    # 坐标轴
    axis.text.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),

    # 标题
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black", margin = margin(b = 20)),

    # 面板
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),

    # 图例
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(4, "cm"),
    legend.margin = margin(l = 20),

    # 边距
    plot.margin = margin(20, 30, 20, 20)
  ) +

  # 标题
  labs(title = "Module-trait relationships")

# 保存ggplot2版本
pdf(file.path("04_Module_Trait_Analysis", "03_Module_Trait_ggplot2_Classic.pdf"),
    width = 8, height = 10)
print(p_classic)
dev.off()

# 输出分析总结
cat("模块-性状关联分析完成！\n")
cat("生成的文件：\n")
cat("1. 01_Module_Trait_Classic_WGCNA.pdf - 经典WGCNA labeledHeatmap风格\n")
cat("2. 02_Module_Trait_Enhanced_Classic.pdf - 增强版（带模块颜色条）\n")
cat("3. 03_Module_Trait_ggplot2_Classic.pdf - ggplot2经典风格\n")

# ------------------------ 基因重要性分析 ---------------------------
update_progress("基因重要性与模块归属度分析")

# 计算基因显著性和模块归属度
geneModuleMembership <- as.data.frame(cor(exprMatrix, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(exprMatrix)))

geneTraitSignificance <- as.data.frame(cor(exprMatrix, clinicalData, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(exprMatrix)))

# 为每个模块和每个分型创建MM vs GS散点图
unique_modules <- unique(moduleColors)

for (module in unique_modules) {
  if (module == "grey") next  # 跳过灰色模块

  # 提取模块基因
  moduleGenes <- (moduleColors == module)
  mm_column <- paste0("ME", module)

  if (mm_column %in% colnames(geneModuleMembership)) {
    # 为每个分型创建散点图
    for (target_trait in unique_clusters) {
      if (target_trait %in% colnames(geneTraitSignificance)) {
        MM <- abs(geneModuleMembership[moduleGenes, mm_column])
        GS <- abs(geneTraitSignificance[moduleGenes, target_trait])

        # 计算相关性
        cor_test <- cor.test(MM, GS, method = "pearson")
        cor_value <- round(cor_test$estimate, 3)
        p_value <- format(cor_test$p.value, scientific = TRUE, digits = 2)

        # 创建数据框
        plot_data <- data.frame(MM = MM, GS = GS)

        # 散点图
        p_scatter <- ggplot(plot_data, aes(x = MM, y = GS)) +
          geom_point(color = module, alpha = 0.7, size = 2) +
          geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
          labs(x = paste("Module Membership in", module, "module"),
               y = paste("Gene Significance for", target_trait),
               title = paste0("Module: ", module, " | Cluster: ", target_trait,
                             "\nCorrelation = ", cor_value, ", P-value = ", p_value)) +
          theme_sci() +
          annotate("text", x = Inf, y = Inf,
                   label = paste("n =", sum(moduleGenes)),
                   hjust = 1.1, vjust = 1.1, size = 4)

        pdf(file.path("05_Gene_Analysis", paste0("MM_vs_GS_", module, "_", target_trait, ".pdf")),
            width = 6, height = 6)
        print(p_scatter)
        dev.off()
      }
    }
  }
}

# ------------------------ 模块可视化 ---------------------------
update_progress("模块可视化与网络图")

# 模块基因数统计
module_sizes <- table(moduleColors)
size_data <- data.frame(
  Module = names(module_sizes),
  GeneCount = as.numeric(module_sizes)
)
size_data$Module <- factor(size_data$Module, levels = size_data$Module[order(size_data$GeneCount, decreasing = TRUE)])

# 模块大小条形图
p_sizes <- ggplot(size_data, aes(x = Module, y = GeneCount, fill = Module)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  scale_fill_identity() +
  labs(title = "Gene Counts per Module",
       x = "Module", y = "Number of Genes") +
  theme_sci() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  geom_text(aes(label = GeneCount), vjust = -0.3, size = 3)

pdf(file.path("07_Module_Visualization", "01_Module_Gene_Counts.pdf"),
    width = figure_width, height = figure_height - 2)
print(p_sizes)
dev.off()

# 模块特征基因表达热图
cat("开始生成模块表达热图...\n")

# 检查是否存在必要的变量
if (!exists("geneModuleMembership") || !exists("moduleColors") || !exists("group_info")) {
  cat("警告：缺少必要的变量，跳过模块热图生成\n")
} else {
  # 过滤掉grey模块，选择前6个最大的模块
  valid_modules <- size_data$Module[size_data$Module != "grey"]
  top_modules <- head(valid_modules, 6)

  cat("将为以下模块生成热图:", paste(top_modules, collapse = ", "), "\n")

  for (module in top_modules) {
    tryCatch({
      # 获取当前模块的基因
      module_genes <- colnames(exprMatrix)[moduleColors == module]
      if (length(module_genes) < 5) {
        cat("模块", module, "基因数少于5，跳过\n")
        next
      }

      cat("处理模块:", module, "，基因数:", length(module_genes), "\n")

      # 选择前50个最相关的基因
      mm_col <- paste0("ME", module)
      if (mm_col %in% colnames(geneModuleMembership)) {
        # 获取当前模块基因的模块归属度
        module_gene_indices <- which(moduleColors == module)
        gene_mm <- abs(geneModuleMembership[module_gene_indices, mm_col])
        names(gene_mm) <- colnames(exprMatrix)[module_gene_indices]

        # 选择前50个最相关的基因
        top_genes <- names(sort(gene_mm, decreasing = TRUE))[1:min(50, length(gene_mm))]

        # 提取表达数据
        module_expr <- t(exprMatrix[, top_genes, drop = FALSE])

        # 检查样本名称匹配
        expr_samples <- colnames(module_expr)
        matched_groups <- group_info[expr_samples]

        # 创建样本注释
        sample_annotation <- data.frame(
          Cluster = factor(matched_groups, levels = unique_clusters)
        )
        rownames(sample_annotation) <- expr_samples

        # 检查是否有NA值
        if (any(is.na(matched_groups))) {
          cat("警告：模块", module, "存在未匹配的样本\n")
        }

        # 注释颜色（使用分型颜色）
        cluster_colors_list <- cluster_colors_map[unique_clusters]
        names(cluster_colors_list) <- unique_clusters

        ann_colors <- list(
          Cluster = cluster_colors_list
        )

        # 生成热图
        pdf_path <- file.path("07_Module_Visualization",
                             paste0("Module_", module, "_Expression_Heatmap.pdf"))

        pdf(pdf_path, width = figure_width + 2, height = figure_width)

        pheatmap(module_expr,
                 annotation_col = sample_annotation,
                 annotation_colors = ann_colors,
                 scale = "row",
                 clustering_distance_rows = "correlation",
                 clustering_distance_cols = "euclidean",
                 show_rownames = FALSE,
                 show_colnames = TRUE,
                 fontsize = 10,
                 fontsize_col = 8,
                 color = colorRampPalette(c("#2166AC", "white", "#D73027"))(100),
                 main = paste("Expression Heatmap for", module, "Module"),
                 border_color = NA)

        dev.off()
        cat("✓ 模块", module, "热图已生成:", pdf_path, "\n")

      } else {
        cat("警告：未找到模块", module, "的特征基因列，跳过\n")
      }

    }, error = function(e) {
      cat("错误：生成模块", module, "热图时出错:", e$message, "\n")
      # 确保PDF设备关闭
      if (dev.cur() > 1) dev.off()
    })
  }
}

# ------------------------ 数据输出 ---------------------------
update_progress("数据输出与文件整理")

# 基因信息汇总
gene_info <- data.frame(
  Gene = colnames(exprMatrix),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

# 添加模块归属度信息
for (mod in colnames(MEs)) {
  gene_info[, paste0("MM_", mod)] <- geneModuleMembership[, mod]
  gene_info[, paste0("MM_pvalue_", mod)] <- MMPvalue[, mod]
}

# 添加基因显著性信息
for (trait in colnames(clinicalData)) {
  gene_info[, paste0("GS_", trait)] <- geneTraitSignificance[, trait]
  gene_info[, paste0("GS_pvalue_", trait)] <- GSPvalue[, trait]
}

# 保存基因信息
write.table(gene_info,
            file = file.path("06_Data_Output", "01_Gene_Module_Information.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 保存模块特征基因
write.table(MEs,
            file = file.path("06_Data_Output", "02_Module_Eigengenes.txt"),
            sep = "\t", row.names = TRUE, quote = FALSE)

# 保存相关性矩阵
write.table(moduleTraitCor,
            file = file.path("06_Data_Output", "03_Module_Trait_Correlations.txt"),
            sep = "\t", row.names = TRUE, quote = FALSE)

write.table(moduleTraitPvalue,
            file = file.path("06_Data_Output", "04_Module_Trait_Pvalues.txt"),
            sep = "\t", row.names = TRUE, quote = FALSE)

# 为每个模块输出基因列表
for (module in unique(moduleColors)) {
  module_genes <- gene_info$Gene[gene_info$Module == module]
  write.table(module_genes,
              file = file.path("06_Data_Output", paste0("05_Module_", module, "_Genes.txt")),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# ------------------------ 分析摘要报告 ---------------------------
update_progress("生成分析摘要报告")

# 创建分析摘要
summary_stats <- list(
  total_genes = ncol(exprMatrix),
  total_samples = nrow(exprMatrix),
  num_clusters = num_clusters,
  cluster_names = paste(unique_clusters, collapse = ", "),
  cluster_sizes = sapply(unique_clusters, function(x) sum(clinicalData[[x]])),
  total_modules = length(unique(moduleColors)),
  soft_threshold = optimalPower,
  largest_module = names(sort(table(moduleColors), decreasing = TRUE))[1],
  largest_module_size = max(table(moduleColors))
)

# 构建分型样本数信息
cluster_info <- paste(sapply(unique_clusters, function(x) {
  paste(x, ":", sum(clinicalData[[x]]), "samples")
}), collapse = "\n")

summary_text <- paste(
  "=== WGCNA分析摘要 ===",
  paste("分析时间:", Sys.time()),
  paste("总基因数:", summary_stats$total_genes),
  paste("总样本数:", summary_stats$total_samples),
  paste("分型数量:", summary_stats$num_clusters),
  paste("分型名称:", summary_stats$cluster_names),
  "各分型样本数:",
  cluster_info,
  paste("检测到的模块数:", summary_stats$total_modules),
  paste("使用的软阈值:", summary_stats$soft_threshold),
  paste("最大模块:", summary_stats$largest_module),
  paste("最大模块基因数:", summary_stats$largest_module_size),
  "",
  "=== 输出文件说明 ===",
  "01_Quality_Control/: 质量控制图表",
  "02_Network_Analysis/: 网络分析图表",
  "03_Module_Detection/: 模块检测图表",
  "04_Module_Trait_Analysis/: 模块-性状关联分析",
  "05_Gene_Analysis/: 基因水平分析（每个模块与每个分型的关联）",
  "06_Data_Output/: 数据输出文件",
  "07_Module_Visualization/: 模块可视化",
  "",
  sep = "\n"
)

writeLines(summary_text, file.path("06_Data_Output", "00_Analysis_Summary.txt"))

# ------------------------ 完成 ---------------------------
end_time <- Sys.time()
total_time <- round(difftime(end_time, start_time, units = "mins"), 2)

cat("\n==== WGCNA分析完成 ====\n")
cat("总用时:", total_time, "分钟\n")
cat("所有结果已保存到相应文件夹中\n")
cat("请查看 06_Data_Output/00_Analysis_Summary.txt 获取详细摘要\n")

# 清理临时变量
rm(list = setdiff(ls(), c("exprMatrix", "MEs", "moduleColors", "clinicalData",
                          "unique_clusters", "group_info", "cluster_colors_map")))
gc()