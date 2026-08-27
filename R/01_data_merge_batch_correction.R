# ==============================================================================
# 01_data_merge_batch_correction.R
# Merge multiple expression matrices and remove batch effects using limma
# Output: Boxplots + advanced PCA plots (with ellipses) + per-sample PCA scores.
# All outputs are written to a timestamped sub-directory.
# ==============================================================================

library(patchwork)
library(data.table)
library(limma)
library(ggplot2)
library(reshape2)
library(tools)
library(ggpubr)
library(RColorBrewer)

# 1. Working directory and output folder
#    All input CSVs are read from the current working directory, and all
#    outputs are written to a timestamped sub-folder inside it.
#    No absolute paths are used so the script is portable across systems.
data_path <- getwd()
out_dirname <- paste0("BatchEffect_output_", format(Sys.time(), "%Y%m%d_%H%M%S"))
output_dir <- file.path(data_path, out_dirname)
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("All outputs will be saved to: ", output_dir, "\n")

# 2. 读取&合并表达数据（CSV格式）
file_list <- list.files(data_path, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("未找到任何csv文件！")
cat("读取到以下文件:\n")
print(file_list)
data_list <- lapply(file_list, function(f){
  dt <- fread(f)
  setnames(dt, 1, "geneSymbol")
  return(dt)
})
merged_data <- Reduce(function(x, y) merge(x, y, by="geneSymbol"), data_list)
cat("合并数据后维度:", dim(merged_data), "\n")
write.csv(
  merged_data, file = file.path(output_dir, "merged_before_batch_removal.csv"),
  row.names = FALSE
)

# 3. 批次名=文件名
batch_names <- basename(file_list)
batch_names <- file_path_sans_ext(batch_names)
sample_counts <- sapply(data_list, function(x) ncol(x) - 1)
batch_info <- unlist(mapply(rep, batch_names, sample_counts))

# 4. 表达矩阵预处理
expr_matrix <- merged_data[, -1, with=FALSE]
expr_matrix <- as.data.frame(expr_matrix)
rownames(expr_matrix) <- merged_data$geneSymbol
expr_matrix <- as.matrix(expr_matrix)
mode(expr_matrix) <- "numeric"
if (length(batch_info) != ncol(expr_matrix)) stop("批次向量与表达矩阵列数不符！")

# 5. 去批次效应
cat("进行批次校正...\n")
corrected_expr <- removeBatchEffect(expr_matrix, batch = batch_info)
corrected_data <- data.frame(geneSymbol = rownames(expr_matrix), corrected_expr, check.names = FALSE)
write.csv(
  corrected_data, file = file.path(output_dir, "merged_after_batch_removal.csv"),
  row.names = FALSE
)

# 6. 简单风格箱线图（PDF保存）
expr_df <- as.data.frame(expr_matrix)
expr_df$gene <- rownames(expr_matrix)
expr_melt <- reshape2::melt(expr_df, id.vars = "gene", variable.name = "Sample", value.name = "Expression")
expr_melt$Batch <- batch_info[match(expr_melt$Sample, colnames(expr_matrix))]

# Nature风格配色
n_batch <- length(unique(batch_info))
nature_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
if(n_batch > length(nature_colors)) {
  nature_colors <- colorRampPalette(nature_colors)(n_batch)
}

p_before <- ggplot(expr_melt, aes(x = Sample, y = Expression, fill = Batch)) +
  geom_boxplot(outlier.size = 0.3, color = NA, lwd = 0) +
  scale_fill_manual(values = nature_colors[1:n_batch]) +
  labs(title = "Before Batch Effect Removal", x = "Sample", y = "Expression") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )
ggsave(file.path(output_dir, "Boxplot_Before_Batch_Removal.pdf"), p_before, width=10, height=5)

corr_expr_df <- as.data.frame(corrected_expr)
corr_expr_df$gene <- rownames(expr_matrix)
corr_expr_melt <- reshape2::melt(corr_expr_df, id.vars = "gene", variable.name = "Sample", value.name = "Expression")
corr_expr_melt$Batch <- batch_info[match(corr_expr_melt$Sample, colnames(expr_matrix))]
p_after <- ggplot(corr_expr_melt, aes(x = Sample, y = Expression, fill = Batch)) +
  geom_boxplot(outlier.size = 0.3, color = NA, lwd = 0) +
  scale_fill_manual(values = nature_colors[1:n_batch]) +
  labs(title = "After Batch Effect Removal", x = "Sample", y = "Expression") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )
ggsave(file.path(output_dir, "Boxplot_After_Batch_Removal.pdf"), p_after, width=10, height=5)

# 7. 高级PCA函数
pca_advanced_plot <- function(mat, batch_info, title, file=NULL){
  data_pca <- prcomp(t(mat), scale. = TRUE)
  pc_var   <- data_pca$sdev^2 / sum(data_pca$sdev^2)  # 方差解释比例
  pca_df <- data.frame(
    Sample = rownames(data_pca$x),
    PC1   = data_pca$x[,1],
    PC2   = data_pca$x[,2],
    Batch = factor(batch_info, levels=unique(batch_info))
  )
  n_batch <- length(unique(batch_info))
  mycol <- brewer.pal(min(n_batch, 8), "Set1")
  if(n_batch > 8) mycol <- colorRampPalette(brewer.pal(9,"Set1"))(n_batch)
  pc1_lab <- paste0("PC1 (", round(pc_var[1]*100, 1),"%)")
  pc2_lab <- paste0("PC2 (", round(pc_var[2]*100, 1),"%)")
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size=4, alpha=0.90) +
    stat_ellipse(aes(fill=Batch), geom="polygon", alpha=0.14, linetype=2, show.legend=FALSE) +
    scale_color_manual(values = mycol) +
    scale_fill_manual(values = mycol) +
    labs(title=title, x=pc1_lab, y=pc2_lab, color="Batch") +
    theme_bw(base_size=14) +
    theme(
      panel.grid.major = element_line(colour="gray90", linetype=2),
      legend.title = element_text(face="bold"),
      legend.background = element_rect(colour="black", fill=NA, size=0.25),
      plot.title = element_text(size=16, face="bold"),
      axis.title = element_text(face="bold")
    )
  if(!is.null(file)){
    ggsave(file, p, width=8, height=6)
  }
  return(p)
}

# 输出PCA PDF
pca_advanced_plot(expr_matrix, batch_info, "PCA Before Batch Effect Removal",
                  file.path(output_dir, "PCA_Advanced_Before_Batch_Removal.pdf"))
pca_advanced_plot(corrected_expr, batch_info, "PCA After Batch Effect Removal",
                  file.path(output_dir, "PCA_Advanced_After_Batch_Removal.pdf"))

# 8. 样本批次表
sample_batch_df <- data.frame(
  Sample = colnames(expr_matrix),
  Batch  = batch_info,
  stringsAsFactors = FALSE
)
write.csv(
  sample_batch_df, file = file.path(output_dir, "sample_batch_info.csv"),
  row.names = FALSE
)

# 9. PCA 得分表（合并前/后）
data_pca_before <- prcomp(t(expr_matrix), scale. = TRUE)
pca_score_before <- data.frame(
  Sample = colnames(expr_matrix),
  PC1 = data_pca_before$x[,1],
  PC2 = data_pca_before$x[,2],
  Batch = batch_info
)
write.csv(
  pca_score_before, file = file.path(output_dir, "PCA_scores_before.csv"),
  row.names = FALSE
)

data_pca_after <- prcomp(t(corrected_expr), scale. = TRUE)
pca_score_after <- data.frame(
  Sample = colnames(expr_matrix),
  PC1 = data_pca_after$x[,1],
  PC2 = data_pca_after$x[,2],
  Batch = batch_info
)
write.csv(
  pca_score_after, file = file.path(output_dir, "PCA_scores_after.csv"),
  row.names = FALSE
)

# 组合箱线图、PCA图
pca_pre  <- pca_advanced_plot(expr_matrix, batch_info, "PCA Before Batch Effect Removal")
pca_post <- pca_advanced_plot(corrected_expr, batch_info, "PCA After Batch Effect Removal")

combined_plot <- (p_before + p_after) / (pca_pre + pca_post) +
  plot_annotation(tag_levels = 'A')  # 自动A,B,C,D

ggsave(file.path(output_dir, "Combined_BatchEffect_Figure.pdf"), combined_plot, width=16, height=10)

cat("流程全部结束，所有文件见: ", output_dir, "\n")
