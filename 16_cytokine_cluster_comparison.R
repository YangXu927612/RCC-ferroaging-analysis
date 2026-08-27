# ==============================================================================
# 16_cytokine_cluster_comparison.R
# Differential expression analysis of inflammatory cytokines between consensus
# clustering subtypes.
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
cytokine_file   <- "Inflammatory Cytokines gene.txt"
cluster_file    <- "geneCluster.txt"

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

# 读取炎症细胞因子基因列表
cytokineGenes <- readLines(cytokine_file)
cytokineGenes <- trimws(cytokineGenes)
cytokineGenes <- cytokineGenes[cytokineGenes != ""]

# 读取聚类分型
cluster_info <- read.table(cluster_file, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
colnames(cluster_info) <- "Cluster"
cluster_info$Cluster <- paste0("C", gsub("C", "", cluster_info$Cluster))

# 匹配炎症因子基因
cytoMatched <- cytokineGenes[cytokineGenes %in% rownames(numExpr)]
cat("炎症因子基因：", length(cytokineGenes), "个，匹配到：", length(cytoMatched), "\n")
if (length(cytoMatched) == 0) stop("炎症因子基因在表达矩阵中均未找到！")

# 匹配样本（聚类文件中的样本）
common_samples <- intersect(colnames(numExpr), rownames(cluster_info))
cat("匹配到的样本数：", length(common_samples), "\n")
if (length(common_samples) == 0) stop("聚类文件与表达矩阵无共同样本！")

# 分型信息
group_vec <- cluster_info[common_samples, "Cluster"]
group_levels <- sort(unique(group_vec))
n_groups <- length(group_levels)
cat("分型数量：", n_groups, "，分别为：", paste(group_levels, collapse=", "), "\n")
cat("各分型样本数：\n")
print(table(group_vec))

# 配色（支持2-4组）
group_colors <- c("#4393C3", "#D6604D", "#66C2A5", "#FC8D62")
group_colors <- group_colors[1:n_groups]
names(group_colors) <- group_levels

# ============================================================
# 第二部分：炎症因子在各分型中的差异箱线图
# ============================================================
cat("\n==== 第二部分：炎症因子差异箱线图 ====\n")

# 提取炎症因子表达（样本×基因）
cytoExpr <- t(numExpr[cytoMatched, common_samples, drop=FALSE])
cytoExpr_df <- as.data.frame(cytoExpr)
cytoExpr_df$Cluster <- group_vec
cytoExpr_df$Sample <- rownames(cytoExpr_df)

# 转长格式
data_long <- melt(cytoExpr_df, id.vars=c("Sample", "Cluster"),
                  variable.name="Cytokine", value.name="Expression")

# 箱线图
boxplot_cyto <- ggboxplot(
  data_long, x="Cytokine", y="Expression", fill="Cluster",
  palette=group_colors,
  xlab="", ylab="Expression Level", legend.title="Cluster", notch=FALSE, width=0.8
) +
  stat_compare_means(
    aes(group=Cluster), label="p.signif",
    symnum.args=list(cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","ns")),
    label.y.npc="top", vjust=-0.5
  ) +
  theme_classic(base_size=14) +
  theme(legend.position="top", axis.text.x=element_text(angle=45, hjust=1),
        axis.line=element_line(color="black", size=1),
        axis.ticks=element_line(color="black", size=1),
        plot.margin=margin(20,10,10,10)) +
  labs(title="Inflammatory Cytokines Expression by Cluster",
       subtitle=paste(sapply(group_levels, function(g) paste0(g, " n=", sum(group_vec==g))), collapse=", ")) +
  coord_cartesian(clip="off")

ggsave("Cytokine_Cluster_Boxplot.pdf", boxplot_cyto, width=max(10, length(cytoMatched)*0.6), height=6)
cat(">> 箱线图已保存为：Cytokine_Cluster_Boxplot.pdf\n")

# p值表
if (n_groups == 2) {
  pvalue_table <- data_long %>%
    group_by(Cytokine) %>%
    do(tidy(wilcox.test(Expression ~ Cluster, data=.))) %>%
    dplyr::select(Cytokine, p.value)
} else {
  pvalue_table <- data_long %>%
    group_by(Cytokine) %>%
    do(tidy(kruskal.test(Expression ~ Cluster, data=.))) %>%
    dplyr::select(Cytokine, p.value)
}
write.csv(pvalue_table, file="Cytokine_Cluster_PValues.csv", row.names=FALSE)
cat(">> p值表已保存\n")

# ============================================================
# 第三部分：炎症因子堆叠条形图（按分型分组）
# ============================================================
cat("\n==== 第三部分：炎症因子堆叠条形图 ====\n")

library(RColorBrewer)

# 样本按分型排序
sample_order <- c()
for (grp in group_levels) {
  sample_order <- c(sample_order, rownames(cytoExpr_df)[cytoExpr_df$Cluster == grp])
}
data_long$Sample <- factor(data_long$Sample, levels=sample_order)

# 配色
cyto_types <- unique(data_long$Cytokine)
nColors <- length(cyto_types)
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)

barplot_cyto <- ggplot(data_long, aes(x=Sample, y=Expression, fill=Cytokine)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=myColors) +
  scale_x_discrete(drop=FALSE) +
  theme_minimal(base_size=18) +
  theme(text=element_text(face="bold"), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major.x=element_blank()) +
  labs(x=NULL, y="Expression Level", fill="Cytokine",
       title="Inflammatory Cytokines Expression Distribution") +
  coord_cartesian(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.05)))

# 动态添加各分型的标注线和标签
seg_colors <- c("#D65DB1", "#0089BA", "#2E8B57", "#FF8C00")
x_start <- 0.5
for (i in seq_along(group_levels)) {
  grp <- group_levels[i]
  grp_count <- sum(cytoExpr_df$Cluster == grp)
  x_end <- x_start + grp_count
  x_mid <- (x_start + x_end) / 2
  barplot_cyto <- barplot_cyto +
    annotate("segment", x=x_start, xend=x_end, y=-0.04, yend=-0.04, color=seg_colors[i], size=5) +
    annotate("text", x=x_mid, y=-0.08, label=grp, color="black", size=7, fontface="bold")
  x_start <- x_end + 1
}

ggsave("Cytokine_barplot_distribution.pdf", barplot_cyto, width=12, height=7.5)
cat(">> 堆叠条形图已保存为：Cytokine_barplot_distribution.pdf\n")

cat("\n==== 全流程执行完毕！====\n")
