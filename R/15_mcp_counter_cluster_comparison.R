# ==============================================================================
# 15_mcp_counter_cluster_comparison.R
# Differential abundance analysis of MCP-counter immune cell populations across
# consensus clustering subtypes.
# ==============================================================================

# ============================ Load required packages =============================
library(MCPcounter)
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
work_dir       <- getwd()                                  # working directory (relative; no absolute paths used)
input_file     <- "merged_after_batch_removal.csv"
cluster_file   <- "geneCluster.txt"

setwd(work_dir)

# ============================================================
# 第一部分：MCP-counter 免疫细胞丰度估算
# ============================================================
cat("==== 第一部分：MCP-counter 免疫细胞丰度估算 ====\n")

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

# 运行MCP-counter免疫细胞丰度估算
cat("正在运行MCP-counter算法...\n")

immune_abundance <- MCPcounter.estimate(
  expression = numExpr,
  featuresType = "HUGO_symbols",
  genes = read.table("genes.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE),
  probesets = read.table("probesets.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
)

# 转换为数据框
immune_df <- as.data.frame(t(immune_abundance))

# 输出结果
write.csv(immune_df, file="MCP-counter-Results.csv", row.names=TRUE)
cat(">> MCP-counter分析完成，结果已保存至 MCP-counter-Results.csv\n")
cat(">> 免疫细胞类型数：", ncol(immune_df), "\n")

# ============================================================
# 第二部分：读取聚类分型信息
# ============================================================
cat("\n==== 第二部分：读取聚类分型 ====\n")

cluster_info <- read.table(cluster_file, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
colnames(cluster_info) <- "Cluster"
cluster_info$Cluster <- paste0("C", gsub("C", "", cluster_info$Cluster))

# 匹配样本
common_samples <- intersect(rownames(immune_df), rownames(cluster_info))
cat("匹配到的样本数：", length(common_samples), "\n")

if (length(common_samples) == 0) stop("聚类文件与MCP-counter结果无共同样本！")

rt <- immune_df[common_samples, , drop=FALSE]
rt$Group <- cluster_info[common_samples, "Cluster"]
rt$Sample <- rownames(rt)

group_levels <- sort(unique(rt$Group))
n_groups <- length(group_levels)
cat("分型数量：", n_groups, "，分别为：", paste(group_levels, collapse=", "), "\n")
cat("各分型样本数：\n")
print(table(rt$Group))

# 配色（支持2-4组）
group_colors <- c("#4393C3", "#D6604D", "#66C2A5", "#FC8D62")
names(group_colors) <- group_levels[1:min(n_groups, 4)]

# ============================================================
# 第三部分：免疫细胞丰度可视化（箱线图 + 堆叠条形图）
# ============================================================
cat("\n==== 第三部分：免疫细胞丰度可视化 ====\n")

# 转为长格式
data_long <- melt(rt, id.vars=c("Sample","Group"), variable.name="Immune", value.name="Abundance")

# 样本按分型排序
sample_order <- c()
for (grp in group_levels) {
  sample_order <- c(sample_order, rownames(rt)[rt$Group == grp])
}
data_long$Sample <- factor(data_long$Sample, levels=sample_order)

# 各组样本数
group_counts <- table(rt$Group)
subtitle_text <- paste(sapply(group_levels, function(g) paste0(g, " n=", group_counts[g])), collapse=", ")

# ---- 箱线图：各分型免疫细胞比较 ----
boxplot_immune <- ggboxplot(
  data_long, x="Immune", y="Abundance", fill="Group",
  palette=group_colors[1:n_groups],
  xlab="", ylab="Abundance Score", legend.title="Cluster", notch=FALSE, width=0.8
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
  labs(title="Immune Cell Abundance by Cluster",
       subtitle=subtitle_text) +
  coord_cartesian(clip="off")

ggsave("immune_abundance_boxplot.pdf", boxplot_immune, width=max(8, n_groups*2.5), height=6)
cat(">> 箱线图已保存为：immune_abundance_boxplot.pdf\n")

# 导出统计表
summary_table <- data_long %>%
  filter(!is.na(Group)) %>%
  group_by(Immune, Group) %>%
  summarise(MeanAbundance=mean(Abundance, na.rm=TRUE), MedianAbundance=median(Abundance, na.rm=TRUE),
            SD=sd(Abundance, na.rm=TRUE), Count=n(), .groups="drop")
write.csv(summary_table, file="boxplot_summary_table.csv", row.names=FALSE)

# p值表（2组用wilcox，多组用kruskal）
if (n_groups == 2) {
  pvalue_table <- data_long %>%
    filter(!is.na(Group)) %>%
    group_by(Immune) %>%
    do(tidy(wilcox.test(Abundance ~ Group, data=.))) %>%
    dplyr::select(Immune, p.value)
} else {
  pvalue_table <- data_long %>%
    filter(!is.na(Group)) %>%
    group_by(Immune) %>%
    do(tidy(kruskal.test(Abundance ~ Group, data=.))) %>%
    dplyr::select(Immune, p.value)
}
write.csv(pvalue_table, file="immune_pvalues.csv", row.names=FALSE)
cat(">> 统计汇总表和p值表已保存\n")

# ---- 堆叠条形图 ----
immune_types <- unique(data_long$Immune)
nColors <- length(immune_types)
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)

barplot_immune <- ggplot(data_long, aes(x=Sample, y=Abundance, fill=Immune)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=myColors) +
  scale_x_discrete(drop=FALSE) +
  theme_minimal(base_size=18) +
  theme(text=element_text(face="bold"), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major.x=element_blank()) +
  labs(x=NULL, y="Relative Abundance", fill="Immune\nCell Type",
       title="Immune Cell Abundance Distribution") +
  coord_cartesian(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.05)))

# 动态添加各分型的标注线和标签
seg_colors <- c("#D65DB1", "#0089BA", "#2E8B57", "#FF8C00")
x_start <- 0.5
for (i in seq_along(group_levels)) {
  grp <- group_levels[i]
  grp_count <- sum(rt$Group == grp)
  x_end <- x_start + grp_count
  x_mid <- (x_start + x_end) / 2
  barplot_immune <- barplot_immune +
    annotate("segment", x=x_start, xend=x_end, y=-0.04, yend=-0.04, color=seg_colors[i], size=5) +
    annotate("text", x=x_mid, y=-0.08, label=grp, color="black", size=7, fontface="bold")
  x_start <- x_end + 1  # +1 留间隔
}

ggsave("barplot_abundance_distribution.pdf", barplot_immune, width=12, height=7.5)
cat(">> 堆叠条形图已保存为：barplot_abundance_distribution.pdf\n")

cat("\n==== 全流程执行完毕！====\n")
