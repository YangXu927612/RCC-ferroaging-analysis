# ==============================================================================
# 12_consensus_clustering.R
# Consensus clustering of RCC samples based on hub candidate gene expression.
# Identifies robust molecular subtypes and outputs cluster labels.
# ==============================================================================

##########################################
#        0. Parameters and thresholds      #
##########################################
# File and directory parameters (relative paths only; no absolute paths used)
gene_exp_file    <- "merged_after_batch_removal.csv"  # input expression matrix (containing both normal and disease samples)
gene_list_file   <- "Intersect_Genes.csv"             # target gene list file
work_dir         <- getwd()                                                                      # working directory (relative path)
output_filename  <- "geneCluster.txt"                  # cluster label output filename
disease_suffix   <- "_tre"                             # suffix used to mark disease samples

# 数据预处理参数
min_row_mean     <- 0                      # 基因表达均值过滤的下限（排除表达量均值 <= 此值的基因）

# 聚类算法相关参数
maximum_k        <- 11                     # 最大聚类类别数
final_cluster_num<- 2                      # 输出分型数（即选第几个K结果）
total_iterations <- 50                     # 重复聚类次数
item_sampling    <- 0.8                    # 每次聚类抽取样本比例
feature_sampling <- 1                      # 每次聚类抽取特征比例
clust_algorithm  <- "km"                   # 聚类算法: "km"=KMeans, "hc"=hierarchical
distance_metric  <- "euclidean"            # 距离度量方法
consensus_seed   <- 2026                 # 共识聚类随机种子
set_rngvalue     <- 2026                 # 全局RNG种子

# 聚类标签映射
labelmap <- c("C1","C2","C3","C4","C5","C6") # 最多支持8类输出标签


##########################################
#        1. 包依赖检测与加载              #
##########################################
cat("步骤1/8: 检查与加载依赖包...\n")
show_progress <- function(step, total=8) {
  pct <- round(step/total*100)
  cat(paste0("进度：", paste(rep("=", step), collapse=""), paste(rep(" ", total-step), collapse=""), " ", pct, "%\n"))
}

if (!requireNamespace("limma", quietly=TRUE)) {
  install.packages("limma")
}
suppressMessages(library(limma))

if (!requireNamespace("ConsensusClusterPlus", quietly=TRUE)) {
  if(!requireNamespace("BiocManager", quietly=TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("ConsensusClusterPlus")
}
suppressMessages(library(ConsensusClusterPlus))
show_progress(1)
library(ggplot2)
library(ggsci)
library(ggsignif)  # 用于显著性星号标注
##########################################
#      2. 工作目录与文件参数检查          #
##########################################
cat("步骤2/8: 检查设置目录和数据文件...\n")

if (!dir.exists(work_dir)) {
  stop("无法找到指定工作目录: ", work_dir)
}
setwd(work_dir)

if (!file.exists(gene_exp_file)) {
  stop("缺少输入文件: ", gene_exp_file)
}
if (!file.exists(gene_list_file)) {
  stop("缺少基因列表文件: ", gene_list_file)
}
show_progress(2)

##########################################
#          3. 数据加载与基础整理           #
##########################################
cat("步骤3/8: 读取与初步检查数据...\n")

raw_data <- tryCatch({
  read.table(gene_exp_file, header=TRUE, sep=",", check.names=FALSE)
}, error=function(e){
  stop("读取基因表达文件失败！")
})

if (ncol(raw_data) < 2) {
  stop("数据表缺少必要的基因/样品列，文件可能损坏！")
}

# 读取目标基因列表
gene_list <- read.csv(gene_list_file, header=TRUE, stringsAsFactors=FALSE)
target_genes <- gene_list[,1]
cat(sprintf("目标基因列表共 %d 个基因。\n", length(target_genes)))

# 设置行名为基因名
raw_matrix <- as.matrix(raw_data)
rownames(raw_matrix) <- as.character(raw_matrix[,1])
expr_data  <- raw_matrix[,-1, drop=FALSE]

# 提取疾病组样本（后缀为 _tre）
tre_cols <- grep(paste0(disease_suffix, "$"), colnames(expr_data), value=TRUE)
if(length(tre_cols) == 0) stop("未找到疾病组样本（后缀为 ", disease_suffix, "）！")
cat(sprintf("疾病组样本数: %d (共 %d 个样本)\n", length(tre_cols), ncol(expr_data)))
expr_data <- expr_data[, tre_cols, drop=FALSE]

# 筛选目标基因（取交集）
matched_genes <- intersect(target_genes, rownames(expr_data))
if(length(matched_genes) == 0) stop("目标基因与表达矩阵无交集，请检查基因名是否一致！")
cat(sprintf("匹配到 %d / %d 个目标基因。\n", length(matched_genes), length(target_genes)))
expr_data <- expr_data[matched_genes, , drop=FALSE]

sample_names <- colnames(expr_data)
gene_names   <- rownames(expr_data)

# 强制为数值型
temp_matrix <- suppressWarnings(matrix(as.numeric(expr_data), nrow=nrow(expr_data)))
dimnames(temp_matrix) <- list(gene_names, sample_names)

if(any(is.na(temp_matrix))){
  warning("表达矩阵中含有NA数据！")
}
if(nrow(temp_matrix)==0 || ncol(temp_matrix)==0){
  stop("表达数据转换后维度为0，请检查表达文件内容！")
}
show_progress(3)

##########################################
#       4. 表达量过滤与合并重复探针        #
##########################################
cat("步骤4/8: 表达矩阵过滤及重复探针合并...\n")
# 表达量过滤
filtered_matrix <- temp_matrix[rowMeans(temp_matrix) > min_row_mean, , drop=FALSE]
if (nrow(filtered_matrix) == 0) stop("根据min_row_mean筛选后无有效基因！")
# 合并重复probe
final_data <- avereps(filtered_matrix)
if(nrow(final_data)==0) stop("avereps处理后为空，原始数据可能异常！")
show_progress(4)

##########################################
#          5. 共识聚类分析                 #
##########################################
cat("步骤5/8: 共识聚类计算...\n")
set.seed(set_rngvalue)
if(maximum_k > ncol(final_data)){
  warning("maximum_k参数大于样本数，注意聚类数量设置。")
}
if(total_iterations < 20){
  warning("total_iterations很小可能导致聚类不稳定。")
}

cluster_result <- tryCatch({
  ConsensusClusterPlus(final_data,
                       maxK        = maximum_k,
                       reps        = total_iterations,
                       pItem       = item_sampling,
                       pFeature    = feature_sampling,
                       title       = work_dir,
                       clusterAlg  = clust_algorithm,
                       distance    = distance_metric,
                       seed        = consensus_seed,
                       plot        = "png")
}, error = function(e) {
  stop("聚类分析出错: ", e$message)
})
show_progress(5)

##########################################
#      6. 聚类结果处理及输出准备          #
##########################################
cat("步骤6/8: 获得聚类结果并重新编码标签...\n")
if(final_cluster_num < 2 | final_cluster_num > maximum_k){
  stop("final_cluster_num设置超限，需在2与最大聚类数之间！")
}

cls_vector <- cluster_result[[final_cluster_num]][["consensusClass"]]
if(is.null(cls_vector)){
  stop("无法检索到指定聚类的聚类标签，请检查前期参数设置！")
}
cls_df <- data.frame(ClusterType=cls_vector)
uniq_levels <- levels(factor(cls_df$ClusterType))
cls_df$ClusterType <- labelmap[match(cls_df$ClusterType, uniq_levels)]

out_df <- rbind(ID=colnames(cls_df), cls_df)
show_progress(6)

##########################################
#                7. 输出                  #
##########################################
cat("步骤7/8: 写出聚类分型标签表...\n")
tryCatch({
  write.table(out_df, file=output_filename, sep="\t", quote=FALSE, col.names=FALSE)
}, error = function(e){
  stop("分型结果写入失败!", e$message)
})
show_progress(7)

# 假设已有下面对象（与上面代码一致）：
# cls_df       # 样本分型标签data.frame，行为样本名
# final_data   # 基因表达矩阵（行为基因，列为样本） 
# cluster_result   # 共识聚类原始结果list
# final_cluster_num   # 当前选择的分型k
# maximum_k          # 最大聚类k

########## 1. 分型统计表 ###########
cat("输出分型统计表...\n")
cluster_stats <- data.frame(
  Cluster = sort(unique(cls_df$ClusterType)),
  Count = as.integer(table(cls_df$ClusterType))
)
cluster_stats$Proportion <- round(cluster_stats$Count / sum(cluster_stats$Count), 4)
write.table(cluster_stats, file = "Cluster_Size_Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########## 2. 当前K的共识矩阵 ###########
cat("输出当前K的共识矩阵...\n")
# 共识矩阵行为样本列为样本，数值0~1表示相似概率
consensus_mat <- cluster_result[[final_cluster_num]]$consensusMatrix
write.table(round(consensus_mat, 4), file = paste0("Consensus_Matrix_K",final_cluster_num,".txt"), sep = "\t", quote = FALSE)

########## 3. 各分型均值表达表 ###########
cat("输出各分型均值表达表 (Feature x Cluster)...\n")
all_clusters <- sort(unique(cls_df$ClusterType))
mean_exp <- sapply(all_clusters, function(grp) {
  idx <- rownames(cls_df)[cls_df$ClusterType == grp]
  rowMeans(final_data[, idx, drop = FALSE])
})
mean_exp <- data.frame(Gene = rownames(final_data), mean_exp)
write.table(mean_exp, file = "Cluster_MeanExp.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########## 4. 不同K的全标签表 ###########
cat("输出全部K分型的标签表...\n")
all_labels <- data.frame(Sample = colnames(final_data))
for (k in 2:maximum_k) {
  class_k <- cluster_result[[k]]$consensusClass
  uniq_lab <- levels(factor(class_k))
  # 用labelmap自动转字母
  all_labels[[paste0("ClusterK",k)]] <- labelmap[match(class_k, uniq_lab)]
}
write.table(all_labels, file = "All_K_Labels.txt", sep = "\t", row.names = FALSE, quote = FALSE)

########## 5. 分型-临床信息合并表（如有临床表） ###########
# 假定有一个clinical.csv ，行为样本名
if(file.exists("clinical.csv")) {
  cat("输出分型-临床合并表...\n")
  clin <- read.csv("clinical.csv", row.names = 1)
  merge_out <- merge(cls_df, clin, by = "row.names", all.x = TRUE)
  write.table(merge_out, file = "Sample_Cluster_Clinical.txt", sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("Step: PCA principal component analysis & advanced visualization...\n")

# 1. Check and align data
if (!("ClusterType" %in% colnames(cls_df))) stop("'ClusterType' column not found in cls_df.")
if (!all(colnames(final_data) %in% rownames(cls_df))) stop("Sample names mismatch between final_data and cls_df.")

# Ensure matching order
cluster_labels <- as.character(cls_df[colnames(final_data), "ClusterType"])

# 2. Transpose for PCA
data4pca <- t(final_data)
pca <- prcomp(data4pca, scale. = TRUE)

# 3. Assemble PCA dataframe & variance explained
expl_var <- round(100 * summary(pca)$importance[2, 1:2], 1) # PC1/2解释方差百分比
pca_df <- as.data.frame(pca$x[, 1:2])
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Subtype <- cluster_labels
pca_df$Sample <- rownames(pca_df)

# 4. Advanced plot using ggplot2
library(ggplot2)
library(ggsci)  # for nice color palettes

set.seed(202406)

# Pick smart palette (如亚型不超9，用nejm/brewer/set1/sci，否则随机)
n_clu <- length(unique(pca_df$Subtype))
if (n_clu <= 9) {
  palette_use <- pal_nejm("default")(n_clu)
} else {
  palette_use <- scales::hue_pal()(n_clu)
}
names(palette_use) <- sort(unique(pca_df$Subtype))

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype)) +
  # 高级点：半透明、边框&大号
  geom_point(size = 4, alpha = 0.78, shape=21, stroke=0.7, aes(fill=Subtype), color="black") +
  scale_color_manual(values=palette_use, name = "Subtype") + 
  scale_fill_manual(values=palette_use, name = "Subtype") +
  theme_bw(base_size = 18) +
  labs(
    title = "PCA of Samples by Subtype",
    x = paste0("PC1 (", expl_var[1], "%)"),
    y = paste0("PC2 (", expl_var[2], "%)")
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face="bold",size=20,hjust=0.5, color="gray25"),
    axis.text = element_text(color="gray30"),
    legend.title=element_text(size=18, face="bold"),
    legend.text=element_text(size=16),
    panel.border=element_rect(color="gray60", size=1.1)
  ) +
  # 椭圆增强分组分布感
  stat_ellipse(aes(group=Subtype, color=Subtype), type="norm", linetype="dashed", size=1.1, alpha=0.5, show.legend=FALSE) +
  guides(color=guide_legend(override.aes=list(size=6)))

# 5. 输出并展示高清图片
ggsave("Cluster_PCA_advanced.png", p, width=7.5, height=6.2, dpi=330)
print(p)
cat("High-level PCA plot saved as Cluster_PCA_advanced.png\n")

# 6. 导出PCA坐标英文标签表
write.table(pca_df, file="PCA_scores_and_subtype.txt", sep="\t", quote=F, row.names=F)
cat("PCA+Subtype scores table output: PCA_scores_and_subtype.txt\n")


cat("全部表格已导出。\n")

# ========================== #
#    筛选差异基因 & 箱线图    #
# ========================== #



# 假定有：final_data（基因x样本），cls_df(有ClusterType)
expr <- final_data
group <- factor(cls_df[colnames(expr), "ClusterType"])

# 差异分析（自动适配2组或多组）
levs <- levels(group)
design <- model.matrix(~0+group)
colnames(design) <- levs
fit <- lmFit(expr, design)

if (length(levs) == 2) {
  contrast_cmd <- paste0("diff = ", levs[2], " - ", levs[1])
  contrast_matrix <- makeContrasts(contrasts=contrast_cmd, levels=design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  deg_table <- topTable(fit2, number=Inf, adjust.method="fdr", sort.by = "P")
} else {
  # 多组：构建所有两两比较的contrast，用F检验筛选任意组间有差异的基因
  pair_contrasts <- combn(levs, 2, function(x) paste0(x[2], " - ", x[1]))
  contrast_matrix <- makeContrasts(contrasts=pair_contrasts, levels=design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  deg_table <- topTable(fit2, number=Inf, adjust.method="fdr", sort.by = "F")
}

# 取"有意义"的adj.P.Val基因（如FDR<0.05，且logFC大于0）
sel_deg <- deg_table[deg_table$adj.P.Val < 0.05, ]
# 只保留前30个
if(nrow(sel_deg)>30){
  sel_deg <- sel_deg[1:30, ]
}
top_genes <- rownames(final_data)
cat(sprintf("将展示前%d个有意义差异基因。\n", length(top_genes)))

# 准备长数据
dat_long <- data.frame(
  Gene = factor(rep(top_genes, each=ncol(expr)), levels=top_genes),
  Expression = as.vector(expr[top_genes, ]),
  Cluster = rep(group, times=length(top_genes))
)

# 统计显著性（自动适配：2组用t检验，多组用Kruskal-Wallis）
stat <- aggregate(Expression ~ Gene + Cluster, data=dat_long, mean)
stat_test <- tapply(top_genes, top_genes, function(g) {
  if (length(levs) == 2) {
    x1 <- expr[g, group==levs[1]]
    x2 <- expr[g, group==levs[2]]
    p <- tryCatch(t.test(x1, x2)$p.value, error=function(e) NA)
  } else {
    gene_vals <- as.numeric(expr[g, ])
    p <- tryCatch(kruskal.test(gene_vals ~ group)$p.value, error=function(e) NA)
  }
  symb <- if(is.na(p)) "" else if(p < 0.001) "***" else if(p < 0.01) "**" else if(p < 0.05) "*" else ""
  return(data.frame(Gene=g, p=p, symbol=symb))
})
starDF <- do.call(rbind, stat_test)

# 箱线图
p <- ggplot(dat_long, aes(x=Gene, y=Expression, fill=Cluster)) +
  geom_boxplot(outlier.size=1, position=position_dodge(width=0.8), width=0.7) +
  scale_fill_manual(values=c("#4393C3", "#D6604D","#8DA0CB","#FC8D62","#66C2A5","#E78AC3")) +  # 自动适配多分型配色
  # 画星号
  geom_text(data=starDF, aes(x=Gene, y=max(dat_long$Expression[dat_long$Gene==Gene])*1.08, label=symbol),
            inherit.aes=FALSE, fontface="bold", color="black", size=6, vjust=0) +
  labs(x="Gene", y="Expression", fill="Cluster") +
  theme_bw(base_size = 17) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    axis.ticks.x = element_line(),
    axis.title = element_text(face="bold"),
    legend.title = element_text(face="bold"),
    strip.text = element_text(face="bold"),
    plot.margin = margin(12,15,5,5)
  )

# 根据基因数量自动调整图宽
ggsave("DEG_boxplot_top30.png", p, width = max(12,0.5*length(top_genes)), height = 6.2, dpi=320)
print(p)
cat("已输出论文式箱线图(DEG_boxplot_top30.png)，如有top>30只显示前30。\n")

sig_deg <- deg_table[deg_table$adj.P.Val < 0.05, ]
# 假设你的差异分析结果：deg_table，行名是gene symbol
file_out <- "DEG_all_results_excel.txt"
# 保存时：col.names=NA 会让第一个单元格空白
write.table(deg_table, file=file_out, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
cat(file_out, "输出成功，首格空白与Excel兼容。\n")

cat("Significant DEGs table (adj.P < 0.05) saved as: DEG_significant_results.txt\n")


cat(sprintf("步骤8/8: 分析完成，分型标签保存在 %s/%s\n", work_dir, output_filename))
show_progress(8)
