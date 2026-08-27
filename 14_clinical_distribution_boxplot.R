# ==============================================================================
# 14_clinical_distribution_boxplot.R
# Stacked bar / boxplot of clinical features across consensus clustering subtypes.
# ==============================================================================

library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(readxl)

# ============================ Parameter settings =============================
work_dir       <- getwd()                                # working directory (relative; no absolute paths used)
cluster_file   <- "geneCluster.txt"
cli_file       <- "clinical.xlsx"

setwd(work_dir)

# ============================================================
# 第一部分：读取数据
# ============================================================
cat("==== 读取数据 ====\n")

# 读取分型结果
Cluster <- read.table(cluster_file, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
colnames(Cluster) <- "Cluster"
cat("分型样本数：", nrow(Cluster), "\n")

# 读取临床数据（xlsx格式）
cli <- as.data.frame(read_excel(cli_file))
rownames(cli) <- cli[[1]]
cli <- cli[, -1, drop=FALSE]
cat("临床样本数：", nrow(cli), "，变量：", paste(colnames(cli), collapse=", "), "\n")

# 匹配样本：分型文件行名带_tre后缀，临床文件不带
cluster_ids <- rownames(Cluster)
cluster_base <- gsub("_tre$", "", cluster_ids)

common_base <- intersect(cluster_base, rownames(cli))
cat("共同样本数：", length(common_base), "\n")

keep_idx <- which(cluster_base %in% common_base)
Cluster <- Cluster[keep_idx, , drop=FALSE]
cluster_base <- cluster_base[keep_idx]

cli <- cli[cluster_base, , drop=FALSE]
rownames(cli) <- rownames(Cluster)

# 保存原始数值Age用于箱线图
age_numeric_raw <- as.numeric(cli$Age)

# Age分组
age_median <- median(as.numeric(cli$Age), na.rm=TRUE)
cli$Age <- ifelse(as.numeric(cli$Age) <= age_median,
                  paste0("<=", age_median),
                  paste0(">", age_median))

# 合并
Type <- cbind(Cluster, cli)
Type$Cluster <- paste0("C", gsub("C", "", Type$Cluster))

cat("Cluster分布：\n")
print(table(Type$Cluster))

# ============================================================
# 第二部分：为每个临床变量绘制堆叠条形图
# ============================================================
cat("\n==== 绘制临床变量分布堆叠条形图 ====\n")

# 配色
bar_colors <- colorRampPalette(brewer.pal(8, "Set2"))(20)

cli_vars <- setdiff(colnames(Type), "Cluster")
plot_list <- list()

for (i in seq_along(cli_vars)) {
  var_name <- cli_vars[i]
  df_tmp <- Type[, c("Cluster", var_name)]
  colnames(df_tmp) <- c("Cluster", "Variable")
  df_tmp <- df_tmp[df_tmp$Variable != "unknow", ]

  # 卡方检验
  tab <- table(df_tmp$Cluster, df_tmp$Variable)
  chi_test <- tryCatch(chisq.test(tab), error = function(e) list(p.value = NA))
  pval <- chi_test$p.value
  p_label <- ifelse(is.na(pval), "NA",
             ifelse(pval < 0.001, "p < 0.001",
             ifelse(pval < 0.01, paste0("p = ", formatC(pval, format="f", digits=3)),
             ifelse(pval < 0.05, paste0("p = ", formatC(pval, format="f", digits=3)),
                    paste0("p = ", formatC(pval, format="f", digits=3))))))

  # 计算比例
  prop_tab <- as.data.frame(prop.table(tab, margin=1))
  colnames(prop_tab) <- c("Cluster", "Category", "Proportion")

  n_cat <- length(unique(prop_tab$Category))
  fill_cols <- bar_colors[1:n_cat]
  names(fill_cols) <- levels(factor(prop_tab$Category))

  p <- ggplot(prop_tab, aes(x = Cluster, y = Proportion, fill = Category)) +
    geom_bar(stat = "identity", width = 0.6, color = "white", size = 0.3) +
    scale_fill_manual(values = fill_cols, name = var_name) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.01)) +
    labs(x = "", y = "Proportion", title = var_name, subtitle = p_label) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "red"),
      axis.text.x = element_text(size = 12, color = "black", face = "bold"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.line = element_line(color = "black", size = 0.8),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(10, 10, 10, 10)
    )

  plot_list[[var_name]] <- p

  # 单独保存
  pdf_name <- paste0("Cluster_", var_name, "_barplot.pdf")
  ggsave(pdf_name, p, width = 4, height = 5)
  cat(">> 已保存：", pdf_name, "  ", p_label, "\n")
}

# ============================================================
# 第三部分：合并所有临床变量为一张图
# ============================================================
cat("\n==== 合并输出 ====\n")

combined <- ggarrange(plotlist = plot_list, ncol = length(plot_list), nrow = 1,
                      common.legend = FALSE)

total_width <- 3.5 * length(plot_list)
ggsave("Cluster_Clinical_Combined.pdf", combined, width = total_width, height = 5.5)
cat(">> 合并图已保存为：Cluster_Clinical_Combined.pdf\n")

# ============================================================
# 第四部分：所有临床变量按Cluster分组的箱线图
# ============================================================
cat("\n==== 绘制临床变量箱线图 ====\n")

# 读取原始临床数据（未分组的）
cli_raw <- as.data.frame(read_excel(cli_file))
rownames(cli_raw) <- cli_raw[[1]]
cli_raw <- cli_raw[, -1, drop=FALSE]
cli_raw <- cli_raw[cluster_base, , drop=FALSE]
rownames(cli_raw) <- rownames(Cluster)

# 合并Cluster信息
cli_raw$Cluster <- Type$Cluster
all_var_names <- setdiff(colnames(cli_raw), "Cluster")

bar_colors <- colorRampPalette(brewer.pal(8, "Set2"))(20)
box_plot_list <- list()
bar_plot_list <- list()

for (var_name in all_var_names) {
  vals <- cli_raw[[var_name]]
  is_numeric <- all(!is.na(suppressWarnings(as.numeric(vals))))

  # ---- 箱线图（数值编码：连续变量直接用，分类变量转数值编码） ----
  if (is_numeric) {
    box_vals <- as.numeric(vals)
  } else {
    box_vals <- as.numeric(factor(vals))
  }
  box_tmp <- data.frame(Cluster = cli_raw$Cluster, Value = box_vals)
  box_tmp <- box_tmp[!is.na(box_tmp$Value), ]

  n_groups <- length(unique(box_tmp$Cluster))
  if (n_groups == 2) {
    wt <- tryCatch(wilcox.test(Value ~ Cluster, data = box_tmp), error = function(e) list(p.value = NA))
  } else {
    wt <- tryCatch(kruskal.test(Value ~ Cluster, data = box_tmp), error = function(e) list(p.value = NA))
  }
  pval_box <- wt$p.value
  p_label_box <- ifelse(is.na(pval_box), "NA",
                 ifelse(pval_box < 0.001, "p < 0.001",
                        paste0("p = ", formatC(pval_box, format="f", digits=3))))

  box_palette <- c("#4393C3", "#D6604D", "#66C2A5", "#FC8D62")[1:n_groups]

  p_box <- ggboxplot(
    box_tmp, x = "Cluster", y = "Value", fill = "Cluster",
    palette = box_palette,
    xlab = "", ylab = var_name,
    width = 0.5, notch = FALSE
  ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "red"),
      axis.text.x = element_text(size = 12, color = "black", face = "bold"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.line = element_line(color = "black", size = 0.8),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = var_name, subtitle = p_label_box)

  box_plot_list[[var_name]] <- p_box
  ggsave(paste0("Cluster_", var_name, "_boxplot.pdf"), p_box, width = 4, height = 5)
  cat(">> 已保存：", paste0("Cluster_", var_name, "_boxplot.pdf"), "\n")

  # ---- 堆叠条形图 ----
  if (is_numeric) {
    # 连续变量按中位数分组
    med_val <- median(as.numeric(vals), na.rm = TRUE)
    cat_vals <- ifelse(as.numeric(vals) <= med_val, paste0("<=", med_val), paste0(">", med_val))
  } else {
    cat_vals <- as.character(vals)
  }
  bar_tmp <- data.frame(Cluster = cli_raw$Cluster, Category = cat_vals)
  bar_tmp <- bar_tmp[!is.na(bar_tmp$Category) & bar_tmp$Category != "unknow" & bar_tmp$Category != "", ]

  if (nrow(bar_tmp) > 0 && length(unique(bar_tmp$Category)) >= 2) {
    tab <- table(bar_tmp$Cluster, bar_tmp$Category)
    chi_test <- tryCatch(chisq.test(tab), error = function(e) list(p.value = NA))
    pval <- chi_test$p.value
    p_label <- ifelse(is.na(pval), "NA",
               ifelse(pval < 0.001, "p < 0.001",
                      paste0("p = ", formatC(pval, format="f", digits=3))))

    prop_tab <- as.data.frame(prop.table(tab, margin=1))
    colnames(prop_tab) <- c("Cluster", "Category", "Proportion")

    n_cat <- length(unique(prop_tab$Category))
    fill_cols <- bar_colors[1:n_cat]
    names(fill_cols) <- levels(factor(prop_tab$Category))

    p_bar <- ggplot(prop_tab, aes(x = Cluster, y = Proportion, fill = Category)) +
      geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
      scale_fill_manual(values = fill_cols, name = var_name) +
      scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.01)) +
      labs(x = "", y = "Proportion", title = var_name, subtitle = p_label) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "red"),
        axis.text.x = element_text(size = 12, color = "black", face = "bold"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.line = element_line(color = "black", size = 0.8),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10)
      )

    bar_plot_list[[var_name]] <- p_bar
    ggsave(paste0("Cluster_", var_name, "_barplot.pdf"), p_bar, width = 4, height = 5)
    cat(">> 已保存：", paste0("Cluster_", var_name, "_barplot.pdf"), "\n")
  }
}

# 合并所有箱线图
if (length(box_plot_list) > 0) {
  combined_box <- ggarrange(plotlist = box_plot_list, ncol = length(box_plot_list), nrow = 1)
  total_w <- 3.5 * length(box_plot_list)
  ggsave("Cluster_Clinical_Boxplot_Combined.pdf", combined_box, width = total_w, height = 5.5)
  cat(">> 合并箱线图已保存为：Cluster_Clinical_Boxplot_Combined.pdf\n")
}

# 合并所有条形图
if (length(bar_plot_list) > 0) {
  combined_bar <- ggarrange(plotlist = bar_plot_list, ncol = length(bar_plot_list), nrow = 1)
  total_w <- 3.5 * length(bar_plot_list)
  ggsave("Cluster_Clinical_Barplot_Combined.pdf", combined_bar, width = total_w, height = 5.5)
  cat(">> 合并条形图已保存为：Cluster_Clinical_Barplot_Combined.pdf\n")
}

cat("\n==== 全流程执行完毕！====\n")
