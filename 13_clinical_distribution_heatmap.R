# ==============================================================================
# 13_clinical_distribution_heatmap.R
# Heatmap of clinical feature distributions across consensus clustering subtypes.
# ==============================================================================

library(pheatmap)
library(limma)
library(RColorBrewer)
library(readxl)

# ============================ Parameter settings =============================
work_dir       <- getwd()                                # working directory (relative; no absolute paths used)
input_file     <- "merged_after_batch_removal.csv"
gene_list_file <- "Intersect_Genes.csv"
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
Cluster <- Cluster[order(Cluster$Cluster), , drop=FALSE]
cat("分型样本数：", nrow(Cluster), "\n")

# 读取临床数据（xlsx格式）
cli <- as.data.frame(read_excel(cli_file))
rownames(cli) <- cli[[1]]
cli <- cli[, -1, drop=FALSE]
cat("临床样本数：", nrow(cli), "，变量：", paste(colnames(cli), collapse=", "), "\n")

# 匹配样本：分型文件行名带_tre后缀，临床文件不带
# 创建映射：去掉_tre后缀用于匹配临床
cluster_ids <- rownames(Cluster)
cluster_base <- gsub("_tre$", "", cluster_ids)

# 找到共同样本
common_base <- intersect(cluster_base, rownames(cli))
cat("共同样本数：", length(common_base), "\n")

# 筛选并对齐
keep_idx <- which(cluster_base %in% common_base)
Cluster <- Cluster[keep_idx, , drop=FALSE]
cluster_base <- cluster_base[keep_idx]

# 对齐临床数据（按分型顺序）
cli <- cli[cluster_base, , drop=FALSE]
rownames(cli) <- rownames(Cluster)  # 统一用带_tre的行名

# Age分组：按中位数分为 <=median 和 >median
age_median <- median(as.numeric(cli$Age), na.rm=TRUE)
cli$Age <- ifelse(as.numeric(cli$Age) <= age_median,
                  paste0("<=", age_median),
                  paste0(">", age_median))
cat("年龄中位数：", age_median, "，分为 <=", age_median, " 和 >", age_median, "\n")

# 合并分型和临床
Type <- cbind(Cluster, cli)
Type$Cluster <- paste0("Cluster", Type$Cluster)

# ============================================================
# 第二部分：临床相关性分析（卡方检验）
# ============================================================
cat("\n==== 临床相关性卡方检验 ====\n")

sigVec <- c("Cluster")
for (clinical in colnames(Type[, 2:ncol(Type)])) {
  data_tmp <- Type[c("Cluster", clinical)]
  colnames(data_tmp) <- c("Cluster", "clinical")
  data_tmp <- data_tmp[data_tmp[, "clinical"] != "unknow", ]
  tableStat <- table(data_tmp)
  stat <- chisq.test(tableStat)
  pvalue <- stat$p.value
  Sig <- ifelse(pvalue < 0.001, "***",
         ifelse(pvalue < 0.01, "**",
         ifelse(pvalue < 0.05, "*", "")))
  sigVec <- c(sigVec, paste0(clinical, Sig))
  cat(sprintf("  %s: p=%.4f %s\n", clinical, pvalue, Sig))
}
colnames(Type) <- sigVec

# ============================================================
# 第三部分：读取表达数据并筛选基因
# ============================================================
cat("\n==== 读取表达数据 ====\n")

# 读取表达矩阵
exprRaw <- read.table(input_file, header=TRUE, sep=",", check.names=FALSE)
exprMat <- as.matrix(exprRaw)
rownames(exprMat) <- exprMat[, 1]
exprVals <- exprMat[, -1, drop=FALSE]
numExpr <- matrix(as.numeric(as.matrix(exprVals)), nrow=nrow(exprVals),
                  dimnames=list(rownames(exprVals), colnames(exprVals)))
numExpr <- avereps(numExpr)

# 读取目标基因
gene_list <- read.csv(gene_list_file, header=TRUE, stringsAsFactors=FALSE)
targetGenes <- as.character(gene_list[, 1])
targetGenes <- targetGenes[targetGenes != ""]

matched_genes <- targetGenes[targetGenes %in% rownames(numExpr)]
cat("目标基因：", paste(matched_genes, collapse=", "), "\n")

# 提取基因表达（基因×样本）
data <- numExpr[matched_genes, , drop=FALSE]

# 只保留分型中的样本
data <- data[, rownames(Type), drop=FALSE]

# ============================================================
# 第四部分：绘制热图
# ============================================================
cat("\n==== 绘制热图 ====\n")

# 定义注释颜色（使用Set3调色板扩展）
colorList <- list()
nColors_total <- sum(sapply(colnames(Type), function(x) length(levels(factor(Type[, x])))))
bioCol <- colorRampPalette(brewer.pal(12, "Set3"))(nColors_total)

j <- 0
for (col_name in colnames(Type)) {
  cliLength <- length(levels(factor(Type[, col_name])))
  cliCol <- bioCol[(j + 1):(j + cliLength)]
  j <- j + cliLength
  names(cliCol) <- levels(factor(Type[, col_name]))
  if ("unknow" %in% names(cliCol)) {
    cliCol["unknow"] <- "grey75"
  }
  colorList[[col_name]] <- cliCol
}

# 绘制热图
pdf("heatmap.pdf", width=10, height=6)
pheatmap(data,
         annotation_col = Type,
         annotation_colors = colorList,
         color = colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         scale = "row",
         show_colnames = FALSE,
         show_rownames = TRUE,
         fontsize = 8,
         fontsize_row = 10,
         fontsize_col = 6,
         border_color = NA)
dev.off()
cat(">> 热图已保存为：heatmap.pdf\n")

cat("\n==== 全流程执行完毕！====\n")
