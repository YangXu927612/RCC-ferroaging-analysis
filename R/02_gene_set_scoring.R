# =========================================================
# 02_gene_set_scoring.R
# Differential expression analysis of research-related pathways
# Method: directly compare mean expression of pathway genes between control and disease.
# =========================================================

# ----------------------------- #
# Step 1: Load required packages
# ----------------------------- #
library(ggpubr)       # publication-ready plots
library(reshape2)     # data reshaping
library(ggplot2)      # general plotting
library(RColorBrewer) # color palettes
library(dplyr)        # data manipulation
library(broom)        # tidy statistical results
library(ggsignif)     # significance markers

# ----------------------------- #
# Step 2: File paths and working directory
#        (all paths are RELATIVE; no absolute paths are used)
# ----------------------------- #
expr_csv    <- "merged_after_batch_removal.csv"
geneset_dir <- "gene set"
# Set working directory to current directory; place input files here before running.
wd_target <- getwd()
if (!dir.exists(wd_target)) stop("Working directory does not exist: ", wd_target)
setwd(wd_target)

# ----------------------------- #
# 第3步：读取表达矩阵
# ----------------------------- #
raw_df <- read.table(expr_csv, header = TRUE, sep = ",",
                     check.names = FALSE, stringsAsFactors = FALSE)
rownames(raw_df) <- raw_df[[1]]
expr_df  <- raw_df[, -1, drop = FALSE]
expr_mat <- apply(expr_df, 2, as.numeric)
rownames(expr_mat) <- rownames(raw_df)

cat("表达矩阵维度:", nrow(expr_mat), "基因 x", ncol(expr_mat), "样本\n")

# ----------------------------- #
# 第4步：读取基因集文件
# ----------------------------- #
txt_files <- list.files(geneset_dir, pattern = "\\.txt$",
                        full.names = TRUE, ignore.case = TRUE)
if (length(txt_files) == 0) stop("未找到基因集文件")

cat("找到", length(txt_files), "个基因集文件\n")

# ----------------------------- #
# 第5步：计算每个通路的平均表达评分
# ----------------------------- #
pathway_scores <- data.frame(row.names = colnames(expr_mat))

for (txt_file in txt_files) {
  pathway_name <- gsub("\\.txt$", "", basename(txt_file), ignore.case = TRUE)

  # 读取基因列表
  genes <- read.table(txt_file, header = TRUE, stringsAsFactors = FALSE)
  gene_vector <- as.character(genes[[1]])
  gene_vector <- gene_vector[gene_vector != "" & !is.na(gene_vector)]

  # 匹配表达矩阵中的基因
  matched_genes <- intersect(gene_vector, rownames(expr_mat))
  cat("  -", pathway_name, ":", length(matched_genes), "/",
      length(gene_vector), "基因匹配\n")

  if (length(matched_genes) >= 2) {
    # 取匹配基因的平均表达量作为通路评分
    pathway_scores[[pathway_name]] <- colMeans(expr_mat[matched_genes, ], na.rm = TRUE)
  } else if (length(matched_genes) == 1) {
    pathway_scores[[pathway_name]] <- expr_mat[matched_genes, ]
  } else {
    cat("    警告：无匹配基因，跳过\n")
  }
}

# ----------------------------- #
# 第6步：添加分组信息
# ----------------------------- #
pathway_scores$Group <- ifelse(
  grepl("_con$", rownames(pathway_scores), ignore.case = TRUE), "Control",
  ifelse(grepl("_tre$", rownames(pathway_scores), ignore.case = TRUE), "Disease", NA)
)
pathway_scores$Sample <- rownames(pathway_scores)

# 过滤掉无法识别组别的样本
pathway_scores <- pathway_scores[!is.na(pathway_scores$Group), ]

cat("\nControl样本数:", sum(pathway_scores$Group == "Control"), "\n")
cat("Disease样本数:", sum(pathway_scores$Group == "Disease"), "\n")

# ----------------------------- #
# 第7步：转长格式
# ----------------------------- #
pathway_cols <- setdiff(colnames(pathway_scores), c("Group", "Sample"))
data_long <- reshape2::melt(pathway_scores, id.vars = c("Sample", "Group"),
                  measure.vars = pathway_cols,
                  variable.name = "Pathway",
                  value.name = "Score")
data_long$Score <- as.numeric(data_long$Score)

# 保存评分结果
write.csv(pathway_scores, file = "Pathway_Mean_Expression_Scores.csv", row.names = FALSE)
cat("通路评分已保存至：Pathway_Mean_Expression_Scores.csv\n")

# ----------------------------- #
# 第8步：为每个通路绘制单独的箱线图
# ----------------------------- #
countControl <- sum(pathway_scores$Group == "Control")
countDisease <- sum(pathway_scores$Group == "Disease")
natureColors <- colorRampPalette(brewer.pal(12, "Set3"))(2)
names(natureColors) <- c("Control", "Disease")

for (pw in pathway_cols) {
  pw_data <- data_long[data_long$Pathway == pw, ]
  y_label <- paste0(pw, " Score")

  boxplot_fig <- ggboxplot(
    pw_data,
    x = "Group",
    y = "Score",
    fill = "Group",
    palette = natureColors,
    xlab = "",
    ylab = y_label,
    legend.title = "Group",
    notch = FALSE,
    width = 0.6
  ) +
    geom_signif(
      comparisons = list(c("Control", "Disease")),
      test = "wilcox.test",
      map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
      textsize = 5,
      vjust = -0.2,
      tip_length = 0.02
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 11, color = "black"),
      legend.title = element_text(size = 12, color = "black", face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 12, color = "black", face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
      plot.margin = margin(20, 10, 10, 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = paste0(pw, " Expression Comparison"),
      subtitle = paste0("Control n=", countControl, ", Disease n=", countDisease)
    ) +
    coord_cartesian(clip = "off")

  pdf_name <- paste0("pathway_", pw, "_comparison.pdf")
  ggsave(pdf_name, boxplot_fig, width = 6, height = 5, dpi = 300)
  cat("已保存：", pdf_name, "\n")
}
cat("所有通路箱线图绑制完成\n")

# ----------------------------- #
# 第9步：导出统计结果
# ----------------------------- #
# 统计汇总表
summary_table <- data_long %>%
  group_by(Pathway, Group) %>%
  summarise(
    MeanScore = mean(Score, na.rm = TRUE),
    MedianScore = median(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  )
write.csv(summary_table, file = "pathway_summary_table.csv", row.names = FALSE)

# p值表
pvalue_table <- data_long %>%
  group_by(Pathway) %>%
  do({
    groups <- unique(.$Group)
    if (length(groups) == 2) {
      tryCatch(
        tidy(wilcox.test(Score ~ Group, data = .)),
        error = function(e) data.frame(p.value = NA)
      )
    } else {
      data.frame(p.value = NA)
    }
  }) %>%
  dplyr::select(Pathway, p.value)
write.csv(pvalue_table, file = "pathway_pvalues.csv", row.names = FALSE)

cat("统计汇总表已保存为：pathway_summary_table.csv\n")
cat("p值表已保存为：pathway_pvalues.csv\n")
cat("分析完成！\n")
