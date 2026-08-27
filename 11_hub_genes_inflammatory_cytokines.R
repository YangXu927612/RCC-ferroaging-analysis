# ==============================================================================
# 11_hub_genes_inflammatory_cytokines.R
# Correlation analysis between hub candidate genes and inflammatory cytokines.
# ==============================================================================

# ============================ Load required packages =============================
library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(broom)
library(limma)
library(pheatmap)
library(scales)

# ============================ Parameter settings =============================
work_dir        <- getwd()                                  # working directory (relative; no absolute paths used)
input_file      <- "merged_after_batch_removal.csv"
gene_list_file  <- "Intersect_Genes.csv"
cytokine_file   <- "Inflammatory Cytokines gene.txt"

setwd(work_dir)

# ============================================================
# 第一部分：读取数据
# ============================================================
cat("==== 第一部分：读取数据 ====\n")

# 读取表达数据
exprRaw <- read.table(input_file, header=TRUE, sep=",", check.names=FALSE)
exprMat <- as.matrix(exprRaw)
rownames(exprMat) <- exprMat[, 1]
exprVals <- exprMat[, -1, drop=FALSE]
rowNamesTmp <- rownames(exprVals)
colNamesTmp <- colnames(exprVals)
numExpr <- matrix(as.numeric(as.matrix(exprVals)), nrow=nrow(exprVals),
                  dimnames=list(rowNamesTmp, colNamesTmp))
numExpr <- avereps(numExpr)
cat("表达矩阵维度：", nrow(numExpr), "x", ncol(numExpr), "\n")

# 读取诊断模型基因列表
gene_list <- read.csv(gene_list_file, header=TRUE, stringsAsFactors=FALSE)
hubGenes <- as.character(gene_list[, 1])
hubGenes <- hubGenes[hubGenes != ""]

# 读取炎症细胞因子基因列表
cytokineGenes <- readLines(cytokine_file)
cytokineGenes <- trimws(cytokineGenes)
cytokineGenes <- cytokineGenes[cytokineGenes != ""]

# 匹配
hubMatched <- hubGenes[hubGenes %in% rownames(numExpr)]
cytoMatched <- cytokineGenes[cytokineGenes %in% rownames(numExpr)]
cat("诊断模型基因：", length(hubGenes), "个，匹配到：", length(hubMatched), "\n")
cat("炎症因子基因：", length(cytokineGenes), "个，匹配到：", length(cytoMatched), "\n")

if (length(hubMatched) == 0) stop("诊断模型基因在表达矩阵中均未找到！")
if (length(cytoMatched) == 0) stop("炎症因子基因在表达矩阵中均未找到！")

# 提取表达值（转置为 样本×基因）
hubExpr  <- t(numExpr[hubMatched, , drop=FALSE])
cytoExpr <- t(numExpr[cytoMatched, , drop=FALSE])

# ============================================================
# 第二部分：Spearman相关性分析
# ============================================================
cat("\n==== 第二部分：Spearman相关性分析 ====\n")

cor_mat      <- matrix(NA, nrow=length(hubMatched), ncol=length(cytoMatched),
                       dimnames=list(hubMatched, cytoMatched))
p_val_matrix <- matrix(NA, nrow=length(hubMatched), ncol=length(cytoMatched),
                       dimnames=list(hubMatched, cytoMatched))

for (i in seq_along(hubMatched)) {
  for (j in seq_along(cytoMatched)) {
    test_result <- cor.test(hubExpr[, i], cytoExpr[, j],
                            method="spearman", use="pairwise.complete.obs")
    cor_mat[i, j]      <- test_result$estimate
    p_val_matrix[i, j] <- test_result$p.value
  }
}

cat("相关性矩阵维度：", nrow(cor_mat), "x", ncol(cor_mat), "\n")

# 创建标签矩阵
labels_mat <- matrix("", nrow=nrow(cor_mat), ncol=ncol(cor_mat),
                     dimnames=dimnames(cor_mat))
for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    stars <- ""
    if (!is.na(p_val_matrix[i, j])) {
      if (p_val_matrix[i, j] < 0.001) stars <- "***"
      else if (p_val_matrix[i, j] < 0.01) stars <- "**"
      else if (p_val_matrix[i, j] < 0.05) stars <- "*"
    }
    cor_val <- ifelse(is.na(cor_mat[i, j]), "", round(cor_mat[i, j], 2))
    labels_mat[i, j] <- paste0(cor_val, "\n", stars)
  }
}

# 导出相关性结果表
corrDF <- data.frame()
for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    corrDF <- rbind(corrDF, data.frame(
      Gene = rownames(cor_mat)[i],
      Cytokine = colnames(cor_mat)[j],
      Correlation = cor_mat[i, j],
      P.Value = p_val_matrix[i, j],
      Significance = ifelse(p_val_matrix[i, j] < 0.05, "Yes", "No")
    ))
  }
}
write.csv(corrDF, file="Gene_Cytokine_Correlation_Results.csv", row.names=FALSE)
cat(">> 相关性结果表已保存\n")

# ============================================================
# 第三部分：相关性热图
# ============================================================
cat("\n==== 第三部分：相关性热图 ====\n")

gene_num   <- nrow(cor_mat)
cyto_num   <- ncol(cor_mat)
pdf_width  <- max(8, 4 + 0.4 * cyto_num)
pdf_height <- max(5, 3 + 0.6 * gene_num)
my_col_fun <- colorRampPalette(c("#6b8e23", "white", "#ba6262"))(100)

cor_mat_plot <- as.matrix(as.data.frame(cor_mat))

pdf("Gene_Cytokine_Correlation_Heatmap.pdf", width=pdf_width, height=pdf_height-3)
pheatmap(
  cor_mat_plot, color=my_col_fun, display_numbers=labels_mat,
  number_color="black", cluster_rows=FALSE, cluster_cols=FALSE,
  fontsize_number=7, fontsize_row=12, fontsize_col=9, border_color="grey90",
  main="Hub Genes - Inflammatory Cytokines Correlation Heatmap",
  angle_col=45
)
dev.off()
cat(">> 相关性热图已保存为：Gene_Cytokine_Correlation_Heatmap.pdf\n")

# ============================================================
# 第四部分：方块气泡图
# ============================================================
cat("\n==== 第四部分：方块气泡图 ====\n")

bubble_df <- corrDF[!is.na(corrDF$P.Value) & !is.na(corrDF$Correlation), ]
bubble_df$AbsCorr <- abs(bubble_df$Correlation)

bubble_df$Sig <- ifelse(bubble_df$P.Value < 0.001, "***",
                 ifelse(bubble_df$P.Value < 0.01, "**",
                 ifelse(bubble_df$P.Value < 0.05, "*", "")))

bubble_df$Cytokine <- factor(bubble_df$Cytokine,
                              levels = rev(sort(unique(bubble_df$Cytokine))))
bubble_df$Gene <- factor(bubble_df$Gene, levels = sort(unique(as.character(bubble_df$Gene))))

n_genes <- length(unique(bubble_df$Gene))
n_cyto  <- length(unique(bubble_df$Cytokine))
fig_width  <- max(6, 3 + n_genes * 1.5)
fig_height <- max(8, 2 + n_cyto * 0.42)

pdf("Gene_Cytokine_SquareBubble.pdf", width = fig_width-2, height = fig_height-5)
p_bubble <- ggplot(bubble_df, aes(x = Gene, y = Cytokine)) +
  geom_point(aes(size = AbsCorr, fill = Correlation),
             shape = 22, color = "grey30", stroke = 0.3) +
  scale_fill_viridis_c(option = "viridis",
                       limits = c(-max(bubble_df$AbsCorr), max(bubble_df$AbsCorr)),
                       oob = squish, name = "Correlation",
                       breaks = pretty(range(bubble_df$Correlation), n = 5)) +
  scale_size_continuous(range = c(2, 10), guide = "none") +
  geom_text(data = bubble_df[bubble_df$Sig != "", ],
            aes(label = Sig), color = "white", size = 4, fontface = "bold", vjust = 0.5) +
  theme(
    panel.background = element_rect(fill = "white", color = "black", size = 1.2),
    panel.grid.major = element_line(color = "black", size = 0.6),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.8),
    legend.position = "right",
    legend.key.height = unit(1.8, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )
print(p_bubble)
dev.off()
cat(">> 方块气泡图已保存为：Gene_Cytokine_SquareBubble.pdf\n")

cat("\n==== 全流程执行完毕！====\n")
