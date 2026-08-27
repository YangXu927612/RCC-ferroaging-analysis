# ==============================================================================
# 09_cibersort_immune_infiltration.R
# Immune infiltration analysis pipeline (CIBERSORT + visualization + gene-immune correlation).
#   - Estimates proportions of 22 immune cell types via CIBERSORT.
#   - Computes Spearman correlations between hub genes and immune cell fractions.
# ==============================================================================

# ============================ Load required packages =============================
library(e1071)
library(parallel)
library(preprocessCore)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(broom)
library(ggsci)
library(limma)
library(pheatmap)

# ============================ Parameter settings =============================
work_dir       <- getwd()                               # working directory (relative; no absolute paths used)
input_file     <- "merged_after_batch_removal.csv"
refer_file     <- "refer.txt"
gene_list_file <- "Intersect_Genes.csv"

setwd(work_dir)

# ============================================================
# 第一部分：CIBERSORT 免疫浸润反卷积
# ============================================================
cat("==== 第一部分：CIBERSORT 免疫浸润反卷积 ====\n")

# 核心SVR算法
coreAlgorithm <- function(refer_matrix, target_vector) {
  num_steps <- 3
  runSVM <- function(step_index) {
    nu_param <- c(0.25, 0.5, 0.75)[step_index]
    svm(refer_matrix, target_vector, type="nu-regression", kernel="linear", nu=nu_param, scale=FALSE)
  }
  if (Sys.info()['sysname'] == 'Windows') {
    models_list <- mclapply(1:num_steps, runSVM, mc.cores=1)
  } else {
    models_list <- mclapply(1:num_steps, runSVM, mc.cores=num_steps)
  }
  rmse_vec <- rep(0, num_steps)
  corr_vec <- rep(0, num_steps)
  for (i in 1:num_steps) {
    current_model <- models_list[[i]]
    weight_vals <- t(current_model$coefs) %*% current_model$SV
    weight_vals[which(weight_vals < 0)] <- 0
    norm_weights <- weight_vals / sum(weight_vals)
    weighted_refer <- sweep(refer_matrix, MARGIN=2, norm_weights, '*')
    est_target <- apply(weighted_refer, 1, sum)
    rmse_vec[i] <- sqrt(mean((est_target - target_vector)^2))
    corr_vec[i] <- cor(est_target, target_vector)
  }
  best_idx <- which.min(rmse_vec)
  best_model <- models_list[[best_idx]]
  best_weights <- t(best_model$coefs) %*% best_model$SV
  best_weights[which(best_weights < 0)] <- 0
  final_weights <- best_weights / sum(best_weights)
  list("final_weights" = final_weights, "best_rmse" = rmse_vec[best_idx], "best_corr" = corr_vec[best_idx])
}

# 置换检验函数
doPermutation <- function(num_perm, refer_matrix, mix_matrix) {
  target_list <- as.list(data.matrix(mix_matrix))
  corr_distribution <- numeric()
  perm_progress <- txtProgressBar(min = 0, max = num_perm, style = 3)
  for (perm_counter in 1:num_perm) {
    permuted_target <- as.numeric(target_list[sample(length(target_list), nrow(refer_matrix), replace = TRUE)])
    if (sd(permuted_target) == 0) {
      std_target <- permuted_target
    } else {
      std_target <- (permuted_target - mean(permuted_target)) / sd(permuted_target)
    }
    perm_result <- coreAlgorithm(refer_matrix, std_target)
    corr_distribution <- c(corr_distribution, perm_result$best_corr)
    setTxtProgressBar(perm_progress, perm_counter)
  }
  close(perm_progress)
  list("null_distribution" = sort(corr_distribution))
}

# CIBERSORT主函数
CIBERSORTModified <- function(refer_file, mix_file, num_perm = 0, do_QN = TRUE) {
  if (!file.exists(refer_file)) stop("参考文件未找到！")
  if (!file.exists(mix_file)) stop("混合文件未找到！")
  refer_data <- data.matrix(read.table(refer_file, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  mix_data   <- data.matrix(read.table(mix_file, header=TRUE, sep=",", row.names=1, check.names=FALSE))
  refer_data <- refer_data[order(rownames(refer_data)), ]
  mix_data   <- mix_data[order(rownames(mix_data)), ]
  if (max(mix_data) < 50) mix_data <- 2^mix_data
  if (do_QN) {
    orig_colnames <- colnames(mix_data); orig_rownames <- rownames(mix_data)
    mix_data <- normalize.quantiles(mix_data)
    colnames(mix_data) <- orig_colnames; rownames(mix_data) <- orig_rownames
  }
  common_genes <- intersect(rownames(refer_data), rownames(mix_data))
  if (length(common_genes) == 0) stop("参考数据与混合数据之间没有共同基因！")
  refer_data <- refer_data[common_genes, , drop=FALSE]
  mix_data   <- mix_data[common_genes, , drop=FALSE]
  refer_data <- (refer_data - mean(refer_data)) / sd(as.vector(refer_data))
  null_dist <- NULL
  if (num_perm > 0) {
    null_dist <- doPermutation(num_perm, refer_data, mix_data)$null_distribution
  }
  output_header <- c("Mixture", colnames(refer_data), "P-value", "Correlation", "RMSE")
  output_matrix <- NULL
  total_samples <- ncol(mix_data)
  sample_progress <- txtProgressBar(min=0, max=total_samples, style=3)
  for (sample_idx in 1:total_samples) {
    curr_sample <- mix_data[, sample_idx]
    std_sample <- if (sd(curr_sample)==0) curr_sample else (curr_sample - mean(curr_sample)) / sd(curr_sample)
    alg_result <- coreAlgorithm(refer_data, std_sample)
    p_val <- NA
    if (!is.null(null_dist)) {
      p_val <- 1 - (which.min(abs(null_dist - alg_result$best_corr)) / length(null_dist))
    }
    curr_result <- c(colnames(mix_data)[sample_idx], alg_result$final_weights, p_val, alg_result$best_corr, alg_result$best_rmse)
    output_matrix <- rbind(output_matrix, curr_result)
    setTxtProgressBar(sample_progress, sample_idx)
  }
  close(sample_progress)
  write.table(rbind(output_header, output_matrix), file="CIBERSORT-Results.csv", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)
  result_numeric <- rbind(output_header, output_matrix)
  result_numeric <- result_numeric[,-1]
  result_numeric <- result_numeric[-1,]
  result_numeric <- matrix(as.numeric(unlist(result_numeric)), nrow=nrow(output_matrix))
  rownames(result_numeric) <- colnames(mix_data)
  colnames(result_numeric) <- c(colnames(refer_data), "P-value", "Correlation", "RMSE")
  return(result_numeric)
}

# 运行CIBERSORT
result_table <- CIBERSORTModified(refer_file, input_file, num_perm=1000, do_QN=TRUE)

# 筛选并输出
filtered_table <- result_table[result_table[,"P-value"] < 1, ]
weight_only <- as.matrix(filtered_table[, 1:(ncol(filtered_table)-3)])
final_output <- rbind(id=colnames(weight_only), weight_only)
write.table(final_output, file="CIBERSORT-Results.csv", sep=",", quote=FALSE, col.names=FALSE)
cat(">> CIBERSORT分析完成，结果已保存至 CIBERSORT-Results.csv\n")

# ============================================================
# 第二部分：免疫浸润可视化（箱线图 + 堆叠条形图）
# ============================================================
cat("\n==== 第二部分：免疫浸润可视化 ====\n")

# 读取CIBERSORT结果
rt <- read.table("CIBERSORT-Results.csv", header=TRUE, sep=",", check.names=FALSE, row.names=1)

# 根据样本后缀判断分组
rt$Group <- ifelse(grepl("_con$", rownames(rt), ignore.case=TRUE), "Control",
                   ifelse(grepl("_tre$", rownames(rt), ignore.case=TRUE), "Treat", NA))
rt$Sample <- rownames(rt)

control_samples <- rownames(rt)[rt$Group == "Control"]
treat_samples   <- rownames(rt)[rt$Group == "Treat"]
all_samples_ordered <- c(control_samples, "gap", treat_samples)

# 转为长格式
data_long <- melt(rt, id.vars=c("Sample","Group"), variable.name="Immune", value.name="Fraction")
data_long$Sample <- factor(data_long$Sample, levels=all_samples_ordered)

countControl <- length(unique(data_long$Sample[data_long$Group == "Control"]))
countTreat   <- length(unique(data_long$Sample[data_long$Group == "Treat"]))

# ---- 箱线图：两组免疫细胞比较 ----
boxplot_immune <- ggboxplot(
  data_long, x="Immune", y="Fraction", fill="Group",
  palette=c("Control"="#FFC0CB", "Treat"="#87CEFA"),
  xlab="", ylab="Fraction", legend.title="Group", notch=FALSE, width=0.8
) +
  stat_compare_means(
    aes(group=Group), label="p.signif",
    symnum.args=list(cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","ns")),
    label.y.npc="top", vjust=-0.5
  ) +
  theme_classic(base_size=14) +
  theme(legend.position="top", axis.text.x=element_text(angle=45, hjust=1),
        axis.line=element_line(color="black", size=1),
        axis.ticks=element_line(color="black", size=1),
        plot.margin=margin(20,10,10,10)) +
  labs(title="Immune Cell Comparison",
       subtitle=paste0(" (Control n=", countControl, ", Treat n=", countTreat, ")")) +
  coord_cartesian(clip="off")

ggsave("immune_diff-points_n.pdf", boxplot_immune, width=8, height=6)
cat(">> 箱线图已保存为：immune_diff-points_n.pdf\n")

# 导出统计表
summary_table <- data_long %>%
  filter(!is.na(Group)) %>%
  group_by(Immune, Group) %>%
  summarise(MeanFraction=mean(Fraction, na.rm=TRUE), MedianFraction=median(Fraction, na.rm=TRUE),
            SD=sd(Fraction, na.rm=TRUE), Count=n(), .groups="drop")
write.csv(summary_table, file="boxplot_summary_table.csv", row.names=FALSE)

pvalue_table <- data_long %>%
  filter(!is.na(Group)) %>%
  group_by(Immune) %>%
  do(tidy(wilcox.test(Fraction ~ Group, data=.))) %>%
  dplyr::select(Immune, p.value)
write.csv(pvalue_table, file="immune_pvalues.csv", row.names=FALSE)
cat(">> 统计汇总表和p值表已保存\n")

# ---- 堆叠条形图 ----
immune_types <- unique(data_long$Immune)
nColors <- length(immune_types)
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)

barplot_immune <- ggplot(data_long, aes(x=Sample, y=Fraction, fill=Immune)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=myColors) +
  scale_x_discrete(drop=FALSE) +
  theme_minimal(base_size=18) +
  theme(text=element_text(face="bold"), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major.x=element_blank()) +
  labs(x=NULL, y="Relative Percent", fill="Immune\nCell Type",
       title="Immune Cell Distribution", subtitle="Control vs. Treat") +
  coord_cartesian(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.05)))

control_count <- length(control_samples)
treat_count   <- length(treat_samples)

barplot_immune <- barplot_immune +
  annotate("segment", x=0.5, xend=control_count+0.5, y=-0.04, yend=-0.04, color="#D65DB1", size=5) +
  annotate("text", x=(control_count)/2+0.5, y=-0.08, label="Control", color="#D65DB1", size=7, fontface="bold") +
  annotate("segment", x=control_count+1.5, xend=control_count+treat_count+1.5, y=-0.04, yend=-0.04, color="#0089BA", size=5) +
  annotate("text", x=control_count+(treat_count)/2+1.5, y=-0.08, label="Treat", color="#0089BA", size=7, fontface="bold")

ggsave("barplot_two_lines_set3_gap.pdf", barplot_immune, width=12, height=7.5)
cat(">> 堆叠条形图已保存为：barplot_two_lines_set3_gap.pdf\n")

# ============================================================
# 第三部分：关键基因与免疫浸润细胞的相关性分析
# ============================================================
cat("\n==== 第三部分：基因-免疫浸润细胞相关性分析 ====\n")

# 读取基因表达数据
exprRaw <- read.table(input_file, header=TRUE, sep=",", check.names=FALSE)
exprMat <- as.matrix(exprRaw)
rownames(exprMat) <- exprMat[, 1]
exprVals <- exprMat[, -1, drop=FALSE]
rowNamesTmp <- rownames(exprVals)
colNamesTmp <- colnames(exprVals)
numExpr <- matrix(as.numeric(as.matrix(exprVals)), nrow=nrow(exprVals),
                  dimnames=list(rowNamesTmp, colNamesTmp))
numExpr <- avereps(numExpr)

# 读取目标基因列表（CSV格式，有表头Gene）
gene_list <- read.csv(gene_list_file, header=TRUE, stringsAsFactors=FALSE)
targetGene <- as.character(gene_list[, 1])
targetGene <- targetGene[targetGene != ""]

# 筛选目标基因
matched_genes <- targetGene[targetGene %in% rownames(numExpr)]
cat("目标基因数：", length(targetGene), "，匹配到：", length(matched_genes), "\n")
exprSelected <- t(numExpr[matched_genes, , drop=FALSE])

# 从CIBERSORT结果提取免疫浸润矩阵（宽格式）
immune_cols <- setdiff(colnames(rt), c("Sample", "Group"))
immuneMat <- as.matrix(rt[, immune_cols])
immuneMat <- apply(immuneMat, 2, as.numeric)
rownames(immuneMat) <- rownames(rt)

# 取共同样本
commonSamples <- intersect(rownames(exprSelected), rownames(immuneMat))
cat("共同样本数：", length(commonSamples), "\n")
exprSelected <- exprSelected[commonSamples, , drop=FALSE]
immuneMat    <- immuneMat[commonSamples, , drop=FALSE]

# 移除标准差为0的免疫细胞列
validImmune <- immuneMat[, apply(immuneMat, 2, sd) > 0, drop=FALSE]
cat("有效免疫细胞类型数：", ncol(validImmune), "\n")

# 计算Spearman相关性
cor_mat     <- matrix(NA, nrow=ncol(exprSelected), ncol=ncol(validImmune),
                      dimnames=list(colnames(exprSelected), colnames(validImmune)))
p_val_matrix <- matrix(NA, nrow=ncol(exprSelected), ncol=ncol(validImmune),
                       dimnames=list(colnames(exprSelected), colnames(validImmune)))

for (i in seq_len(ncol(exprSelected))) {
  for (j in seq_len(ncol(validImmune))) {
    test_result <- cor.test(exprSelected[,i], validImmune[,j], method="spearman", use="pairwise.complete.obs")
    cor_mat[i,j]     <- test_result$estimate
    p_val_matrix[i,j] <- test_result$p.value
  }
}

# 创建标签矩阵（相关系数 + 显著性星号）
labels_mat <- matrix("", nrow=nrow(cor_mat), ncol=ncol(cor_mat),
                     dimnames=dimnames(cor_mat))
for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    stars <- ""
    if (!is.na(p_val_matrix[i,j])) {
      if (p_val_matrix[i,j] < 0.001) stars <- "***"
      else if (p_val_matrix[i,j] < 0.01) stars <- "**"
      else if (p_val_matrix[i,j] < 0.05) stars <- "*"
    }
    cor_val <- ifelse(is.na(cor_mat[i,j]), "", round(cor_mat[i,j], 2))
    labels_mat[i,j] <- paste0(cor_val, "\n", stars)
  }
}

# 调试输出
cat("cor_mat 维度：", nrow(cor_mat), "x", ncol(cor_mat), "\n")
cat("cor_mat 行名：", paste(rownames(cor_mat), collapse=", "), "\n")
cat("labels_mat 维度：", nrow(labels_mat), "x", ncol(labels_mat), "\n")

# 绘制相关性热图
gene_num   <- nrow(cor_mat)
immune_num <- ncol(cor_mat)
pdf_width  <- max(8, 4 + 0.5 * immune_num)
pdf_height <- max(5, 3 + 0.6 * gene_num)
my_col_fun <- colorRampPalette(c("#6b8e23", "white", "#ba6262"))(100)

# 确保cor_mat转为data.frame再转回matrix，避免pheatmap的维度问题
cor_mat_plot <- as.matrix(as.data.frame(cor_mat))

pdf("Gene_Immune_Correlation_Heatmap.pdf", width=pdf_width, height=pdf_height)
pheatmap(
  cor_mat_plot, color=my_col_fun, display_numbers=labels_mat,
  number_color="black", cluster_rows=FALSE, cluster_cols=FALSE,
  fontsize_number=7, fontsize_row=12, fontsize_col=9, border_color="grey90",
  main="Hub Genes - Immune Cell Correlation Heatmap",
  angle_col=45
)
dev.off()
cat(">> 相关性热图已保存为：Gene_Immune_Correlation_Heatmap.pdf\n")

# 导出相关性结果表
corrDF <- data.frame()
for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    corrDF <- rbind(corrDF, data.frame(
      Gene=rownames(cor_mat)[i], ImmuneCell=colnames(cor_mat)[j],
      Correlation=cor_mat[i,j], P.Value=p_val_matrix[i,j],
      Significance=ifelse(p_val_matrix[i,j] < 0.05, "Yes", "No")
    ))
  }
}
write.csv(corrDF, file="Gene_Immune_Correlation_Results.csv", row.names=FALSE)
cat(">> 相关性结果表已保存为：Gene_Immune_Correlation_Results.csv\n")

# ============================================================
# ============================================================
# 第四部分：方块气泡图（Square Bubble Chart）基因-免疫细胞相关性
# ============================================================
cat("\n==== 第四部分：方块气泡图 ====\n")

library(ggplot2)
library(scales)

# 准备数据（直接用corrDF，免疫细胞名保持原样）
bubble_df <- corrDF[!is.na(corrDF$P.Value) & !is.na(corrDF$Correlation), ]
bubble_df$AbsCorr <- abs(bubble_df$Correlation)

# 显著性标记
bubble_df$Sig <- ifelse(bubble_df$P.Value < 0.001, "***",
                        ifelse(bubble_df$P.Value < 0.01, "**",
                               ifelse(bubble_df$P.Value < 0.05, "*", "")))

# 调试：打印数据确认
cat("bubble_df 行数：", nrow(bubble_df), "\n")
cat("基因：", paste(unique(bubble_df$Gene), collapse=", "), "\n")
cat("免疫细胞：", paste(unique(bubble_df$Immune), collapse=", "), "\n")

# 免疫细胞按名称排序
bubble_df$Immune <- factor(bubble_df$Immune,
                           levels = rev(sort(unique(bubble_df$Immune))))
bubble_df$Gene <- factor(bubble_df$Gene, levels = sort(unique(as.character(bubble_df$Gene))))

# 图尺寸
n_genes  <- length(unique(bubble_df$Gene))
n_immune <- length(unique(bubble_df$Immune))
fig_width  <- max(6, 3 + n_genes * 1.5)
fig_height <- max(8, 2 + n_immune * 0.42)

pdf("Gene_Immune_SquareBubble.pdf", width = fig_width-2, height = fig_height-4)
p_bubble <- ggplot(bubble_df, aes(x = Gene, y = Immune)) +
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
    plot.margin = ggplot2::margin(10, 10, 10, 10)   # 唯一修改处
  )
print(p_bubble)
dev.off()
cat(">> 方块气泡图已保存为：Gene_Immune_SquareBubble.pdf\n")

cat("\n==== 全流程执行完毕！====\n")