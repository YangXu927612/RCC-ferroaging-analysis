# ==============================================================================
# 03_differential_expression_analysis.R
# LIMMA-based differential expression analysis between RCC and adjacent normal tissues.
# Generates volcano plots, heatmaps, PCA plots, Venn diagrams and intersection
# with ferro-aging-related gene sets.
# ==============================================================================

# Set working directory to current directory (no absolute path).
setwd(getwd())

threshold_logFC <- 0.5       # log2折叠变化阈值
threshold_adjP  <- 0.05    # 调整后P值阈值
max_display_genes <- 50    # 热图中每侧（上调和下调）展示的基因最大数量

# 记录 Step 1 开始时间
step1_start <- Sys.time()
cat("[Step 1] 开始运行...\n")

# 载入必要的 R 包 -----------------------------------------------------
library(limma)    # 用于线性模型和差异表达分析
library(ggplot2)
library(ggrepel)
library(pheatmap) # 用于绘制热图
library(RColorBrewer)
library(grid)
library(ComplexHeatmap)
library(circlize)

# 记录 Step 1 结束时间
step1_end <- Sys.time()
cat("[Step 1] 运行完毕，耗时:", 
    round(difftime(step1_end, step1_start, units="secs"), 2), "秒\n\n")



# ========================== 2. 参数设定及文件路径配置 ==========================
step2_start <- Sys.time()
cat("[Step 2] 开始参数设定和文件路径配置...\n")



# 定义输入文件路径（请根据实际情况修改）---------------------------------------
file_expr  <- "merged_after_batch_removal.csv"  # 基因表达矩阵文件
 

step2_end <- Sys.time()
cat("[Step 2] 运行完毕，耗时:", 
    round(difftime(step2_end, step2_start, units="secs"), 2), "秒\n\n")
# ========================== 读取表达矩阵并保证数值型 ==========================
file_expr <- "merged_after_batch_removal.csv" # 修改成你的实际数据文件名

# 建议自动检测分隔符
tmp_head <- readLines(file_expr, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")

# 读取原始表达文件，行为基因名，列为样本
expr_raw <- read.table(file_expr, header=TRUE, sep=sep, check.names=FALSE, stringsAsFactors=FALSE)
rownames(expr_raw) <- expr_raw[,1]
expr_mat <- expr_raw[,-1, drop=FALSE]

# 强制转为matrix并确保是数值型
expr_mat <- as.matrix(expr_mat)
expr_mat <- apply(expr_mat, 2, as.numeric) # 逐列强制类型
rownames(expr_mat) <- rownames(expr_raw)
colnames(expr_mat) <- colnames(expr_raw)[-1] # 保持列名
combined_expr <- expr_mat

# 检查所有元素皆为数值
if(!is.numeric(expr_mat)) stop("表达矩阵仍存在非数值列！")
cat("表达矩阵维度:", dim(expr_mat)[1], "基因 x", dim(expr_mat)[2], "样本\n")

# ========================== 自动判断分组信息 ==========================
sample_names <- colnames(expr_mat)
group_info <- ifelse(grepl("_con$", sample_names, ignore.case = TRUE), "Control",
                     ifelse(grepl("_tre$", sample_names, ignore.case = TRUE), "Disease", "Unknown"))
if(any(group_info == "Unknown")) stop("分组未知的样本存在，请检查样本名后缀！")
num_ctrl <- sum(group_info == "Control")
num_treat <- sum(group_info == "Disease")
cat("分组情况：Control =", num_ctrl, ", Disease =", num_treat, "\n")

# ========================== 差异表达分析（LIMMA） ==========================
group_labels <- factor(group_info, levels = c("Control", "Disease"))
design_mat <- model.matrix(~0 + group_labels)
colnames(design_mat) <- c("Control", "Disease")
fit <- lmFit(expr_mat, design_mat)
contrast_mat <- makeContrasts(Disease - Control, levels = design_mat)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)
all_diff_results <- topTable(fit2, adjust.method = "fdr", number = Inf)
cat("差异分析完成！发现FDR<0.05基因数：", sum(all_diff_results$adj.P.Val < 0.05), "\n")

# 可选：输出全部基因的差异分析表
write.table(cbind(Gene=rownames(all_diff_results), all_diff_results),
            file="DE_results.csv", sep=",", quote=FALSE, row.names=FALSE)


# 可选：写出到文件
# write.csv(all_diff_results, file="DEG_results.csv")

step5_end <- Sys.time()
#cat("[Step 5] 用时：", round(as.numeric(difftime(step5_end, step5_start, units = "secs")), 2), "秒\n")


# ------------------------- Step 6: 筛选显著差异表达基因并输出包含必要信息的表格 -------------------------
step6_start <- Sys.time()
cat("[Step 6] 开始筛选显著差异表达基因并输出包含必要信息的表格...\n")

# 根据设定的阈值筛选显著上调或下调的基因
significant_DEGs <- all_diff_results[with(all_diff_results,
                                          (abs(logFC) > threshold_logFC & adj.P.Val < threshold_adjP)), ]

# 整理输出格式，添加 "Gene" 列显示基因名称，
# 此外默认的 topTable 输出已经包含以下信息：
#   - logFC: 基因的对数折叠变化
#   - AveExpr: 各样本的平均表达值
#   - t: 经过经验贝叶斯调整后的 t 统计量
#   - P.Value: 原始 p 值
#   - adj.P.Val: 调整后 p 值（FDR校正）
#   - B: 差异表达的对数 odds
output_DEGs <- cbind(Gene = rownames(significant_DEGs), significant_DEGs)

# 计算标准误差 (SE)
# 根据 t 统计量公式： t = logFC / SE  =>  SE = |logFC/t|
# 当 t 不为 0 时计算，否则返回 NA
SE <- ifelse(as.numeric(output_DEGs[, "t"]) != 0,
             abs(as.numeric(output_DEGs[, "logFC"]) / as.numeric(output_DEGs[, "t"])),
             NA)
output_DEGs <- cbind(output_DEGs, SE = SE)

# 如果有额外的注释信息（如 EntrezID、GeneDescription、Chromosome 等），
# 可通过 merge() 将外部注释表与结果合并，示例如下：
# annotation_table <- read.table("gene_annotation.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
# output_DEGs <- merge(output_DEGs, annotation_table, by.x="Gene", by.y="GeneSymbol", all.x=TRUE)

# 可选：调整列顺序，使结果更直观（例如把 Gene, logFC, SE, AveExpr, t, P.Value, adj.P.Val, B 放在前面）
output_DEGs <- as.data.frame(output_DEGs, stringsAsFactors = FALSE)

# 新增Regulation列，显示上下调信息
output_DEGs$Regulation <- ifelse(as.numeric(output_DEGs$logFC) > 0, "Up", "Down")

desired_order <- c("Gene", "logFC", "Regulation", "SE", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
existing_cols <- intersect(desired_order, colnames(output_DEGs))
output_DEGs <- output_DEGs[, c(existing_cols, setdiff(colnames(output_DEGs), existing_cols))]

# 将输出表格写入CSV文件，确保所有信息都展示
write.csv(output_DEGs, file = "DE_significant_genes.csv", row.names = FALSE)

step6_end <- Sys.time()
cat("[Step 6] 运行完毕，耗时:", round(difftime(step6_end, step6_start, units = "secs"), 2), "秒\n\n")


# ========================== 7. 自动识别gene set目录下的txt文件并计算交集 ==========================
step7_start <- Sys.time()
cat("[Step 7] 开始自动识别基因集文件并计算交集...\n")

# 检查并安装ggvenn包
if (!require("ggvenn", quietly = TRUE)) {
  cat("安装ggvenn包...\n")
  install.packages("ggvenn")
  library(ggvenn)
}

# 自动识别 gene set 目录下所有txt文件
gene_set_dir <- "gene set"
gene_set_files <- list.files(gene_set_dir, pattern = "\\.txt$", full.names = TRUE)
cat("在", gene_set_dir, "目录下发现", length(gene_set_files), "个基因集文件\n")

if (length(gene_set_files) == 0) {
  stop("未在 gene set 目录下找到任何txt文件！")
}

# 获取差异基因列表
deg_genes <- rownames(significant_DEGs)

# 存储所有基因集的结果
all_gene_sets <- list()       # 每个基因集的基因列表
all_intersect_results <- list() # 每个基因集与DEG的交集

for (gs_file in gene_set_files) {
  # 获取基因集名称（去掉.txt后缀）
  gs_name <- tools::file_path_sans_ext(basename(gs_file))
  cat("\n--- 处理基因集:", gs_name, "---\n")

  # 读取基因集
  gs_data <- read.table(gs_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

  # 获取基因列（第一列或名为"Genes"/"Gene"的列）
  if ("Genes" %in% colnames(gs_data)) {
    gs_genes <- unique(gs_data$Genes)
  } else if ("Gene" %in% colnames(gs_data)) {
    gs_genes <- unique(gs_data$Gene)
  } else {
    gs_genes <- unique(gs_data[[1]])
  }

  # 移除NA和空值
  gs_genes <- gs_genes[!is.na(gs_genes) & gs_genes != ""]

  # 计算交集
  gs_intersect <- intersect(gs_genes, deg_genes)

  cat("  基因集数量:", length(gs_genes), "\n")
  cat("  差异基因数量:", length(deg_genes), "\n")
  cat("  交集基因数量:", length(gs_intersect), "\n")

  all_gene_sets[[gs_name]] <- gs_genes
  all_intersect_results[[gs_name]] <- gs_intersect
}

# 兼容后续代码：使用最后一个基因集作为默认（用于火山图标记等）
research_genes <- all_gene_sets[[length(all_gene_sets)]]
intersect_genes <- all_intersect_results[[length(all_intersect_results)]]

step7_end <- Sys.time()
cat("\n[Step 7] 运行完毕，耗时:",
    round(difftime(step7_end, step7_start, units="secs"), 2), "秒\n\n")


# ========================== 8. 绘制渐变色风格火山图 ==========================
step8_start <- Sys.time()
cat("[Step 8] 开始绘制渐变色风格火山图（标记交集基因）...\n")

# 准备数据
volcano_data <- all_diff_results
volcano_data$Gene <- rownames(volcano_data)

# 选择要标注的基因：使用交集基因
# 筛选出在交集中的差异基因
intersect_deg_data <- volcano_data[volcano_data$Gene %in% intersect_genes, ]

# 如果交集基因不足40个，则全部展示
if (nrow(intersect_deg_data) <= 40) {
  label_genes_gradient <- intersect_deg_data$Gene
  cat("交集基因不足40个，全部标记:", length(label_genes_gradient), "个基因\n")
} else {
  # 超过40个，则分别展示上调和下调排名前20的基因
  # 上调基因：logFC > 0，按-log10(adj.P.Val)排序
  up_intersect <- intersect_deg_data[intersect_deg_data$logFC > 0, ]
  up_intersect <- up_intersect[order(-(-log10(up_intersect$adj.P.Val))), ]

  # 下调基因：logFC < 0，按-log10(adj.P.Val)排序
  down_intersect <- intersect_deg_data[intersect_deg_data$logFC < 0, ]
  down_intersect <- down_intersect[order(-(-log10(down_intersect$adj.P.Val))), ]

  # 取前20个（如果不足20则全部取）
  up_label <- head(up_intersect$Gene, min(20, nrow(up_intersect)))
  down_label <- head(down_intersect$Gene, min(20, nrow(down_intersect)))

  label_genes_gradient <- unique(c(up_label, down_label))
  cat("交集基因超过40个，标记上调前", length(up_label), "个和下调前", length(down_label), "个，共", length(label_genes_gradient), "个基因\n")
}
volcano_data$label_gene <- ifelse(volcano_data$Gene %in% label_genes_gradient, volcano_data$Gene, "")

# 计算x轴和y轴范围
x_limit <- max(abs(volcano_data$logFC), na.rm = TRUE) * 1.1
y_max <- max(-log10(volcano_data$adj.P.Val), na.rm = TRUE) * 1.05

# 计算-log10(p)的范围，用于设置气泡大小图例的breaks
y_values <- -log10(volcano_data$adj.P.Val)
y_values <- y_values[is.finite(y_values)]
y_range <- range(y_values, na.rm = TRUE)
# 根据数据范围自动生成合适的breaks
size_breaks <- pretty(y_range, n = 5)
size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= max(y_values)]

# 计算上调和下调基因数量
up_count <- sum(volcano_data$logFC > threshold_logFC & volcano_data$adj.P.Val < threshold_adjP)
down_count <- sum(volcano_data$logFC < -threshold_logFC & volcano_data$adj.P.Val < threshold_adjP)

# 绘制渐变色火山图
gradient_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val))) +
  # 绘制点，颜色按logFC渐变，大小按-log10(p值)
  geom_point(aes(color = logFC, size = -log10(adj.P.Val)), alpha = 0.85) +
  # 为交集基因添加黑色外框标记
  geom_point(data = subset(volcano_data, label_gene != ""),
             aes(size = -log10(adj.P.Val)),
             shape = 21, fill = NA, color = "black", stroke = 1.2) +
  # 颜色渐变：深蓝-青绿-黄绿-黄-橙-红
  scale_color_gradientn(
    colors = c("#2B83BA",   # 深蓝
               "#5AAE61",   # 青绿
               "#A6D96A",   # 黄绿
               "#FFFFBF",   # 淡黄
               "#FEE08B",   # 黄
               "#FDAE61",   # 橙
               "#D7191C"),  # 红
    values = scales::rescale(c(-6, -4, -2, 0, 2, 4, 6)),
    limits = c(-x_limit, x_limit),
    name = "log2FC"
  ) +
  # 点大小映射 - 动态设置breaks
  scale_size_continuous(
    range = c(0.8, 7),
    name = "-log10(p_val)",
    breaks = size_breaks,
    labels = as.character(round(size_breaks, 0))
  ) +
  # 添加阈值线
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC),
             linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = -log10(threshold_adjP),
             linetype = "dashed", color = "grey50", linewidth = 0.6) +
  # 添加基因标签 - 只对交集基因子集标注，确保全部显示
  geom_text_repel(
    data = subset(volcano_data, label_gene != ""),
    aes(label = label_gene),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey40",
    segment.size = 0.3,
    fontface = "italic",
    color = "black",
    show.legend = FALSE
  ) +
  # 添加"Down"和"Up"标注箭头（Down箭头朝左，Up箭头朝右）
  annotate("segment", x = -x_limit * 0.55, xend = -x_limit * 0.9,
           y = y_max * 0.97, yend = y_max * 0.97,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           color = "#2B83BA", linewidth = 1.2) +
  annotate("text", x = -x_limit * 0.72, y = y_max * 0.97,
           label = paste0("Down (", down_count, ")"), color = "#2B83BA", fontface = "bold",
           size = 4.5, vjust = -0.8) +
  annotate("segment", x = x_limit * 0.55, xend = x_limit * 0.9,
           y = y_max * 0.97, yend = y_max * 0.97,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           color = "#D7191C", linewidth = 1.2) +
  annotate("text", x = x_limit * 0.72, y = y_max * 0.97,
           label = paste0("Up (", up_count, ")"), color = "#D7191C", fontface = "bold",
           size = 4.5, vjust = -0.8) +
  # 添加p值阈值标注
  annotate("text", x = x_limit * 0.98, y = -log10(threshold_adjP),
           label = paste0("p = ", threshold_adjP),
           hjust = 1, vjust = -0.5, size = 3.5, color = "grey30") +
  # 设置坐标轴
  scale_x_continuous(limits = c(-x_limit, x_limit), expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0.02, 0)) +
  # 标题和标签
  labs(
    title = "Volcano Plot",
    x = "avg_log2FC",
    y = "-log10(p_val)"
  ) +
  # 主题设置
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.2),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.spacing.y = unit(0.3, "cm"),
    plot.margin = margin(15, 15, 10, 10)
  ) +
  # 图例设置
  guides(
    color = guide_colorbar(
      title = "avg_log2FC",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 1.2,
      barheight = 10,
      frame.colour = "black",
      frame.linewidth = 0.5,
      ticks.colour = "black",
      ticks.linewidth = 0.5,
      order = 1
    ),
    size = guide_legend(
      title = "-log10(p_val)",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, color = "grey30"),
      order = 2
    )
  )

# 保存渐变色火山图
ggsave("DE_volcano_gradient.pdf", gradient_volcano, width = 9, height = 7.5, dpi = 300)

step8_end <- Sys.time()
cat("[Step 8] 运行完毕，耗时:", round(difftime(step8_end, step8_start, units="secs"), 2), "秒\n\n")

# ========================== 9. PCA 主成分分析 ==========================
step9_start <- Sys.time()
cat("[Step 9] 开始优化 PCA 主成分分析...\n")

# 确保所需包已加载
# 计算 PCA
pca_result <- prcomp(t(combined_expr), scale. = TRUE)

# 构建 PCA 数据框
pca_df <- data.frame(
  Sample = colnames(combined_expr),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = factor(c(rep("Control", num_ctrl), rep("Disease", num_treat)))  # 分组信息
)

# 计算 PCA 贡献率（用于坐标轴标签）
pca_var <- pca_result$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

# 绘制 PCA 图 - 新风格：方形点、填充椭圆、无标签、右侧图例
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  # 填充椭圆（置信区间）
  stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.2, linewidth = 1) +
  # 方形点
  geom_point(size = 3, shape = 15) +
  # 自定义颜色：Disease蓝色方块，Control橙红色圆点
  scale_color_manual(values = c("Control" = "#E64B35",    # 橙红色
                                "Disease" = "#4DBBD5")) + # 蓝色
  scale_fill_manual(values = c("Control" = "#E64B35",     # 橙红色填充
                               "Disease" = "#4DBBD5")) + # 蓝色填充
  # 修改坐标轴标签
  labs(x = paste0("PCA1(", pca_var_perc[1], "%)"),
       y = paste0("PCA2(", pca_var_perc[2], "%)"),
       color = "Group",
       fill = "Group") +
  # 美化主题
  theme_bw(base_size = 14) +
  theme(
    plot.title    = element_blank(),
    axis.title    = element_text(face = "bold", size = 14, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    legend.title  = element_text(size = 12, face = "bold"),
    legend.text   = element_text(size = 11),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.line     = element_blank()
  )

# 保存 PCA 图
ggsave("DE_PCA.pdf", pca_plot, width = 7, height = 6, dpi = 300)

step9_end <- Sys.time()
cat("[Step 9] PCA 分析优化完成，耗时:", 
    round(difftime(step9_end, step9_start, units="secs"), 2), "秒\n\n")

# 使用真实的差异分析结果替代模拟数据
# 根据logFC和p值排序
gene_data_real <- all_diff_results[order(all_diff_results$logFC, decreasing = TRUE),]

# 创建渐变色
color_palette <- scale_color_gradientn(
  colors = rev(brewer.pal(9, "RdYlBu")),
  values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1))
)

# 绘制基因排序点图
plot <- ggplot(gene_data_real, aes(x = seq_along(logFC), y = logFC)) +
  geom_point(aes(color = adj.P.Val, size = abs(logFC)), alpha = 0.8) +  # 使用p值作为颜色，logFC的绝对值作为大小
  color_palette +  # 使用渐变色来表示p值
  scale_size_continuous(range = c(1, 3)) +  # 调整点的大小范围
  labs(
    title = "Gene Ranking Dotplot", 
    x = "Ranking of Differentially Expressed Genes", 
    y = "Log2 Fold Change",
    color = "Adjusted P-value",
    size = "Log2 Fold Change"
  ) +
  theme_minimal(base_size = 16) +  # 使用最简洁的背景，调整字体大小
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "#2F4F4F", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold", color = "#2F4F4F"),
    axis.text = element_text(size = 14, color = "#555555"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

# 保存为PDF，调整尺寸和分辨率
ggsave("DE_gene_ranking_dotplot.pdf", plot, width = 10, height = 8, dpi = 300)


# 假定你的差异分析结果为 all_diff_results，行名为基因名
de_table <- all_diff_results

# 使用交集基因作为标记基因
label_genes <- intersect_genes

# 添加gene列到全结果表
de_table$Gene <- rownames(de_table)

# 增加一个新列：为交集基因标记gene name，否则空
de_table$label <- ifelse(de_table$Gene %in% label_genes, de_table$Gene, "")

# 颜色分组（可选，可按你的Significance或重新定义）
de_table$Group <- "Not Significant"
de_table$Group[de_table$logFC > threshold_logFC & de_table$adj.P.Val < threshold_adjP] <- "Up regulated"
de_table$Group[de_table$logFC < -threshold_logFC & de_table$adj.P.Val < threshold_adjP] <- "Down regulated"

# 绘制火山图
volcano_plot <- ggplot(de_table, aes(x=logFC, y=-log10(adj.P.Val), color=Group)) +
  geom_point(alpha=0.7, size=2) +
  geom_text_repel(data = subset(de_table, label != ""),
                  aes(label=label),
                  size=3,
                  max.overlaps=Inf,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.color="grey50",
                  color = "black",
                  show.legend = FALSE) +
  scale_color_manual(values = c("Up regulated"="#FF4500",
                                "Down regulated"="#1E90FF",
                                "Not Significant"="#808080")) +
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC), linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(threshold_adjP), linetype="dashed", color="black") +
  labs(title="Volcano Plot with Intersection Genes Labeled",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value") +
  theme_minimal(base_size = 14) +
  theme(plot.title=element_text(face="bold", hjust=0.5, color="#2F4F4F"),
        axis.line = element_line(color = "black", linewidth = 0.5))

# 保存
ggsave("DE_volcano_with_labels.pdf", volcano_plot, width=8, height=8, dpi=300)
# ========================== 全部步骤完成 ==========================

# 假设你的显著性标记列为 Group, 可根据你的实际列名适配
group_counts <- table(de_table$Group)
up_count   <- ifelse("Up regulated"   %in% names(group_counts), group_counts["Up regulated"], 0)
down_count <- ifelse("Down regulated" %in% names(group_counts), group_counts["Down regulated"], 0)
not_count  <- ifelse("Not Significant"%in% names(group_counts), group_counts["Not Significant"], 0)
# 自定义图例名字，把数量拼进去
legend_labels <- c(
  paste0("Up regulated (", up_count, ")"),
  paste0("Down regulated (", down_count, ")"),
  paste0("Not Significant (", not_count, ")")
)
# labels的顺序要和scale_color_manual的values顺序一致
volcano_plot <- ggplot(de_table, aes(x=logFC, y=-log10(adj.P.Val), color=Group)) +
  geom_point(alpha=0.7, size=2) +
  geom_text_repel(data = subset(de_table, label != ""),
                  aes(label=label),
                  size=3,
                  max.overlaps=Inf,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.color="grey50",
                  color = "black",
                  show.legend = FALSE) +
  scale_color_manual(
    values = c("Up regulated" = "#FF4500",
               "Down regulated" = "#1E90FF",
               "Not Significant" = "#808080"),
    breaks = c("Up regulated", "Down regulated", "Not Significant"),
    labels = legend_labels
  ) +
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC), linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(threshold_adjP), linetype="dashed", color="black") +
  labs(title="Volcano Plot with Intersection Genes Labeled",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value") +
  theme_minimal(base_size = 14) +
  theme(plot.title=element_text(face="bold", hjust=0.5, color="#2F4F4F"),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5))

# 保存
ggsave("DE_volcano_with_labels.pdf", volcano_plot, width=8, height=8, dpi=300)

# 差异基因名，以data.frame形式输出TXT，表头为gene
diff_gene_list <- data.frame(gene = rownames(significant_DEGs))
write.table(diff_gene_list,
            file = "DEG_geneList.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

# ========================== 9.5. 为每个基因集输出交集分析结果和Venn图 ==========================
step9_5_start <- Sys.time()
cat("[Step 9.5] 开始为每个基因集输出交集分析结果和Venn图...\n")

for (gs_name in names(all_gene_sets)) {
  cat("\n--- 基因集:", gs_name, "---\n")

  gs_genes <- all_gene_sets[[gs_name]]
  gs_intersect <- all_intersect_results[[gs_name]]

  # 创建Venn图的基因列表
  gene_list_venn <- list()
  gene_list_venn[[gs_name]] <- gs_genes
  gene_list_venn[["DEG"]] <- deg_genes

  # 绘制Venn图
  venn_file <- paste0("DEG_", gs_name, "_Venn.pdf")
  pdf(file = venn_file, width = 8, height = 8)
  print(ggvenn(gene_list_venn,
         show_percentage = TRUE,
         stroke_color = "white",
         stroke_size = 0.5,
         fill_color = c("#1E90FF", "#FF8C00"),
         set_name_color = c("#1E90FF", "#FF8C00"),
         set_name_size = 6,
         text_size = 4.5))
  dev.off()
  cat("  Venn图已保存:", venn_file, "\n")

  # 保存交集基因列表（CSV格式）
  csv_file <- paste0("Intersection_DEG_", gs_name, ".csv")
  write.csv(data.frame(Gene = gs_intersect), file = csv_file, row.names = FALSE)
  cat("  交集基因列表(CSV)已保存:", csv_file, "\n")

  # 保存交集基因列表（TXT格式）
  txt_file <- paste0("Intersection_DEG_", gs_name, ".txt")
  write.table(data.frame(Gene = gs_intersect), file = txt_file,
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  交集基因列表(TXT)已保存:", txt_file, "\n")

  # 输出统计摘要
  cat("  基因集数量:", length(gs_genes), "\n")
  cat("  差异基因数量:", length(deg_genes), "\n")
  cat("  交集基因数量:", length(gs_intersect), "\n")
  if (length(gs_genes) > 0) {
    cat("  交集率(vs基因集):", round(length(gs_intersect)/length(gs_genes)*100, 2), "%\n")
  }
  if (length(deg_genes) > 0) {
    cat("  交集率(vs差异基因):", round(length(gs_intersect)/length(deg_genes)*100, 2), "%\n")
  }
}

step9_5_end <- Sys.time()
cat("[Step 9.5] 运行完毕，耗时:",
    round(difftime(step9_5_end, step9_5_start, units="secs"), 2), "秒\n\n")

# ========================== 10. 绘制带分布密度图的热图 ==========================
step11_start <- Sys.time()
cat("[Step 11] 开始绘制带分布密度图的热图...\n")

# 检查是否有足够的差异基因
if(nrow(significant_DEGs) == 0) {
  cat("警告：未发现显著差异基因，跳过热图绘制\n")
} else {
  # 1. 准备绘图数据
  ordered_DEGs <- significant_DEGs[order(as.numeric(as.vector(significant_DEGs$logFC))), ]
  ordered_gene_names <- rownames(ordered_DEGs)
  total_DEG_count <- length(ordered_gene_names)

  if (total_DEG_count > (max_display_genes * 2)) {
    selected_gene_set <- ordered_gene_names[c(1:max_display_genes,
                                              (total_DEG_count - max_display_genes + 1):total_DEG_count)]
  } else {
    selected_gene_set <- ordered_gene_names
  }
  heatmap_expr <- combined_expr[selected_gene_set, ]

  # 2. 行标准化（Z-score标准化）
  heatmap_expr_scaled <- t(scale(t(heatmap_expr)))

  # 3. 样本注释信息
  sample_annotation <- data.frame(Group = factor(c(rep("Control", num_ctrl),
                                                   rep("Disease", num_treat))))
  rownames(sample_annotation) <- colnames(combined_expr)

  # 4. 分组颜色设置
  annotation_colors <- list(
    Group = c("Control" = "#66C2A5",    # 浅绿
              "Disease" = "#FC8D62")   # 浅橙
  )

  # 5. 创建顶部注释：分布密度图 + 分组注释
  ha_top <- HeatmapAnnotation(
    Distribution = anno_density(heatmap_expr_scaled, type = "heatmap",
                                which = "column", height = unit(2, "cm")),
    Group = sample_annotation$Group,
    col = annotation_colors,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 12)
  )

  # 6. 热图颜色方案（蓝-白-红）
  color_palette <- colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026"))

  # 7. 为基因名添加上下调标记
  # 获取选中基因的logFC值
  gene_logFC <- significant_DEGs[selected_gene_set, "logFC"]
  # 根据logFC添加标记：疾病组上调(logFC>0)加(Up)，下调(logFC<0)加(Down)
  gene_labels <- paste0(selected_gene_set,
                        ifelse(gene_logFC > 0, " (Up)", " (Down)"))
  names(gene_labels) <- selected_gene_set

  # 8. 创建ComplexHeatmap对象
  ht <- Heatmap(
    heatmap_expr_scaled,
    name = "Expression\n(Z-score)",
    col = color_palette,
    top_annotation = ha_top,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    row_labels = gene_labels[rownames(heatmap_expr_scaled)],
    column_title = paste("Differential Expression Heatmap with Distribution\n",
                        "Control: ", num_ctrl, " | Disease: ", num_treat),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    )
  )

  # 8. 保存到PDF
  pdf("DE_heatmap_with_distribution.pdf", width = 12, height = 10)
  tryCatch({
    draw(ht)
    cat("[Step 11] 带分布密度图的热图已保存: DE_heatmap_with_distribution.pdf\n")
  }, error = function(e) {
    cat("绘制ComplexHeatmap时出现错误:", e$message, "\n")
    cat("尝试使用print()方法...\n")
    tryCatch({
      print(ht)
      cat("[Step 11] 带分布密度图的热图已保存: DE_heatmap_with_distribution.pdf\n")
    }, error = function(e2) {
      cat("ComplexHeatmap绘制失败，可能是包未正确安装\n")
    })
  })
  dev.off()
}

step11_end <- Sys.time()
cat("[Step 11] 运行完毕，耗时:",
    round(difftime(step11_end, step11_start, units="secs"), 2), "秒\n\n")



cat("所有步骤已成功完成！\n")


# ========================== 13. 交集基因热图 ==========================
cat("[Step 13] 开始绘制交集基因热图...\n")

# 读取交集基因
intersect_gene_file <- "Intersection_DEG_ferro-aging.txt"
intersect_gene_data <- read.table(intersect_gene_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
plot_genes <- intersect_gene_data[[1]]
plot_genes <- plot_genes[plot_genes %in% rownames(combined_expr)]

if (length(plot_genes) > 0) {
  # 提取交集基因的表达矩阵
  heatmap_mat <- combined_expr[plot_genes, , drop = FALSE]

  # Z-score标准化（按行）
  heatmap_scaled <- t(scale(t(heatmap_mat)))

  # 分组注释
  ha_col <- HeatmapAnnotation(
    Group = factor(c(rep("Control", num_ctrl), rep("Disease", num_treat))),
    col = list(Group = c("Control" = "#2166AC", "Disease" = "#B2182B")),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 11)
  )

  # 热图颜色
  col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

  ht_intersect <- Heatmap(
    heatmap_scaled,
    name = "Z-score",
    col = col_fun,
    top_annotation = ha_col,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 10, fontface = "italic"),
    column_title = "Intersection Genes Heatmap",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )

  pdf("Intersection_genes_heatmap.pdf", width = 10, height = max(4, length(plot_genes) * 0.45))
  draw(ht_intersect)
  dev.off()
  cat("交集基因热图已保存: Intersection_genes_heatmap.pdf\n")
} else {
  cat("警告：交集基因在表达矩阵中未找到，跳过热图\n")
}

# ========================== 14. 交集基因箱线图 ==========================
cat("[Step 14] 开始绘制交集基因箱线图...\n")

if (!require("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr")
  library(ggpubr)
}

if (length(plot_genes) > 0) {
  # 构建长格式数据
  group_vec <- factor(c(rep("Control", num_ctrl), rep("Disease", num_treat)),
                      levels = c("Control", "Disease"))

  box_list <- list()
  for (g in plot_genes) {
    expr_vals <- as.numeric(combined_expr[g, ])
    box_list[[g]] <- data.frame(
      Gene = g,
      Expression = expr_vals,
      Group = group_vec
    )
  }
  box_df <- do.call(rbind, box_list)
  box_df$Gene <- factor(box_df$Gene, levels = plot_genes)

  # 从limma差异分析结果中提取adj.P.Val作为显著性标记（与DE_significant_genes.csv一致）
  p_vals <- sapply(plot_genes, function(g) {
    if (g %in% rownames(all_diff_results)) {
      all_diff_results[g, "adj.P.Val"]
    } else {
      NA
    }
  })

  format_pval <- function(p) {
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("ns")
  }
  p_labels <- sapply(p_vals, format_pval)

  # 计算每个基因的y轴最大值，用于放置显著性标记
  y_positions <- sapply(plot_genes, function(g) {
    max(as.numeric(combined_expr[g, ])) * 1.08
  })

  anno_df <- data.frame(
    Gene = factor(plot_genes, levels = plot_genes),
    label = p_labels,
    y = y_positions
  )

  # 所有基因在同一个x轴上，x轴为Gene，颜色区分Control/Disease
  box_plot <- ggplot(box_df, aes(x = Gene, y = Expression, fill = Group)) +
    geom_boxplot(outlier.size = 0.8, width = 0.65, alpha = 0.85,
                 position = position_dodge(0.75)) +
    scale_fill_manual(values = c("Control" = "#2166AC", "Disease" = "#B2182B")) +
    # 添加显著性标注
    geom_text(data = anno_df, aes(x = Gene, y = y, label = label),
              inherit.aes = FALSE, size = 4, fontface = "bold") +
    labs(x = "", y = "Expression", fill = "Group") +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(size = 11, color = "black", face = "bold.italic",
                                 angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(face = "bold", size = 13),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1)
    )

  plot_width <- max(6, length(plot_genes) * 1.2 + 2)
  ggsave("Intersection_genes_boxplot.pdf", box_plot,
         width = plot_width-5, height = 5.5, dpi = 300)
  cat("交集基因箱线图已保存: Intersection_genes_boxplot.pdf\n")
} else {
  cat("警告：无交集基因可绘制箱线图\n")
}

cat("\n[Step 13-14] 交集基因热图和箱线图绘制完成！\n")

