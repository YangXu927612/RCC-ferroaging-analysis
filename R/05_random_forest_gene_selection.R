# ==============================================================================
# 05_random_forest_gene_selection.R
# Random Forest-based feature importance analysis for gene selection.
# Steps: data loading, preprocessing, model construction, feature selection,
# and visualization.
# ==============================================================================

#----------------------------
# Step 1: Load required R packages
#----------------------------
# Working directory (relative paths only)
workDir <- getwd()

# 自定义参数：基因重要性评分筛选阈值
importance_threshold <- 2   # 输出评分大于此值的基因
library(randomForest)
library(ggplot2)
library(viridis)

# 设置随机种子以确保结果的可重复性
set.seed(2025)  # 修改了种子值

# 初始化进度条
totalSteps <- 10   # 总共10个步骤
pb <- txtProgressBar(min = 0, max = totalSteps, style = 3)
step <- 0  # 当前进度计数器

#----------------------------
# 第2步：设置工作目录和读取数据文件
#----------------------------

step <- step + 1  # 进度更新
setTxtProgressBar(pb, step)  # 更新进度条

# 指定数据文件名称（请确认文件在指定路径下存在）
dataFileName <- "merged_after_batch_removal.csv"

if (dir.exists(workDir)) {               # 判断工作目录是否存在
  setwd(workDir)                         # 如果存在，则设置工作目录
} else {
  dir.create(workDir, recursive = TRUE)   # 如果不存在，则创建目录
  setwd(workDir)
}

# 自动识别 Intersection_ 前缀的txt文件作为基因列表
gene_list_files <- list.files(workDir, pattern = "^Intersection_.*\\.txt$", full.names = TRUE)
if (length(gene_list_files) == 0) {
  stop("未在工作目录下找到 Intersection_ 前缀的txt文件！")
}
geneListFile <- gene_list_files[1]
message(sprintf("自动识别基因列表文件: %s", basename(geneListFile)))

# 读取基因列表（自动判断是否有表头）
geneListRaw <- readLines(geneListFile)
geneListRaw <- trimws(geneListRaw)
geneListRaw <- geneListRaw[geneListRaw != ""]
if (tolower(geneListRaw[1]) %in% c("gene", "genes", "genename", "gene_name")) {
  geneList <- geneListRaw[-1]
} else {
  geneList <- geneListRaw
}
geneList <- unique(geneList)
message(sprintf("共读取 %d 个基因", length(geneList)))

# 判断数据文件是否存在，若不存在则终止执行
if (!file.exists(dataFileName)) {
  stop("数据文件不存在，请检查文件名或路径！")  # 输出错误提示并停止程序
}

# 读取数据文件（假设CSV文件，带标题，第一列为行名）
rawData <- read.table(dataFileName, header = TRUE, sep = ",",
                      check.names = FALSE, row.names = 1)

# 根据基因列表筛选基因
availableGenes <- intersect(geneList, rownames(rawData))
message(sprintf("在表达数据中找到 %d / %d 个基因", length(availableGenes), length(geneList)))
if (length(availableGenes) == 0) {
  stop("指定的基因列表中没有在表达数据中找到的基因！请检查基因名称是否匹配。")
}
rawData <- rawData[availableGenes, , drop = FALSE]

# 复制一份数据用于冗余备份
backupData <- rawData

#----------------------------
# 第3步：数据预处理
#----------------------------

step <- step + 1  # 更新进度条
setTxtProgressBar(pb, step)

# 转置数据，使每行代表一个样本，每列代表一个基因
exprData <- t(rawData)

# 判断转置后数据是否为空
if (nrow(exprData) == 0) {
  stop("转置后的数据为空，请检查输入数据！")
}

# 从样本名称中提取组别信息，假设样本名形如 "sample01_Control" 或 "sample02_Treat"
sampleLabels <- row.names(exprData)
# 使用正则表达式提取下划线后的部分作为组别
groupInfo <- gsub(".*_(.*)$", "\\1", sampleLabels)
# 如果提取结果中存在NA或空字符，则用默认值替代
groupInfo[groupInfo == ""] <- "Unknown"
# 在 Step3 数据预处理后，加上这一行
colnames(exprData) <- make.names(colnames(exprData))

#----------------------------
# 第4步：构建初始随机森林模型
#----------------------------

step <- step + 1  # 更新进度条
setTxtProgressBar(pb, step)

# 构建初始随机森林模型：使用所有基因作为预测变量，groupInfo作为响应变量
# 使用500棵树来构建模型
rfModel <- randomForest(as.factor(groupInfo) ~ ., data = as.data.frame(exprData), ntree = 500)

# 检查模型输出是否正常
if (is.null(rfModel)) {
  stop("随机森林模型构建失败，请检查数据格式！")
}

#----------------------------
# 第5步：保存随机森林模型错误率图
#----------------------------

step <- step + 1  # 更新进度条
setTxtProgressBar(pb, step)

# 绘制并保存模型错误率图至PDF文件
pdf(file = "RF_error_plot.pdf", width = 6, height = 6)  # 开启PDF设备
plot(rfModel, main = "Random Forest Error Rate", lwd = 2, 
     col = c("black", "red", "blue"), type = "l", lty = c(1,2,3))  # 绘制错误率曲线
# 添加图例说明
legend("topright", legend = c("Overall OOB Error", "Class 1 OOB Error", "Class 2 OOB Error"),
       col = c("black", "red", "blue"), lwd = 2, lty = c(1,2,3), cex = 0.8, box.lwd = 0.5)
dev.off()  # 关闭PDF设备，保存图像

#----------------------------
# 第6步：确定最佳决策树数量
#----------------------------

step <- step + 1  # 更新进度条
setTxtProgressBar(pb, step)

# 从模型中找出错误率最低时的树数量
bestTreeCount <- which.min(rfModel$err.rate[, 1])
# 如果最佳树数量小于某个阈值，则提示并调整
if (bestTreeCount < 50) {
  bestTreeCount <- 50  # 强制最少使用50棵树
}

# 输出最佳树数量（打印到控制台）
print(paste("最佳决策树数量为:", bestTreeCount))

#----------------------------
# 第7步：重新训练随机森林模型
#----------------------------

step <- step + 1  # 更新进度条
setTxtProgressBar(pb, step)

# 使用最佳树数量重新训练随机森林模型
rfOptimized <- randomForest(as.factor(groupInfo) ~ ., data = as.data.frame(exprData), ntree = bestTreeCount)

# 判断重新训练模型是否成功
if (is.null(rfOptimized)) {
  stop("重新训练的随机森林模型构建失败！")
}

##----------------------------
# 第8步：提取和筛选基因重要性
#----------------------------

step <- step + 1  # 更新进度条
setTxtProgressBar(pb, step)

# 获取模型中每个变量的重要性评分
rawImportance <- importance(rfOptimized)

# 按照Gini指数降序排列基因
if (!is.null(dim(rawImportance))) {
  # 如果返回的是矩阵，则提取 "MeanDecreaseGini" 列
  giniValues <- rawImportance[,"MeanDecreaseGini"]
} else {
  # 如果返回的是向量，则直接使用
  giniValues <- rawImportance
}
# 对基因按照重要性进行排序
sortedGini <- sort(giniValues, decreasing = TRUE)

# 选择前10个最重要的基因
topGenes <- names(sortedGini)[1:min(10, length(sortedGini))]

# 如果基因数量低于10个，改为按评分大于阈值进行筛选
if (length(sortedGini) < 10) {
  topGenes <- names(sortedGini[sortedGini > importance_threshold])
  message(sprintf("基因数量不足10个，改为按评分>%g筛选，共 %d 个基因", importance_threshold, length(topGenes)))
  if (length(topGenes) == 0) {
    message("警告：没有评分大于阈值的基因，保留评分最高的基因")
    topGenes <- names(sortedGini)
  }
}

# 按自定义阈值筛选基因并输出
filteredGenes <- names(sortedGini[sortedGini > importance_threshold])
message(sprintf("评分大于 %g 的基因共 %d 个", importance_threshold, length(filteredGenes)))

if (length(filteredGenes) == 0) {
  message("警告：没有评分大于阈值的基因，将保留评分最高的基因")
  filteredGenes <- names(sortedGini)[1]
}

# 将筛选出的特征基因保存到文本文件中
write.table(filteredGenes, file = "Top_RF_Genes.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
message(sprintf("已输出 %d 个基因到 Top_RF_Genes.txt（阈值: %g）", length(filteredGenes), importance_threshold))

#----------------------------
# 第9步：绘制基因重要性气泡图（使用ggplot2和viridis）
#----------------------------

step <- step + 1  # 更新进度条
setTxtProgressBar(pb, step)

# 构造数据框用于绘图：这里提取所有基因及其重要性
impDF <- data.frame(Gene = names(sortedGini), 
                    Importance = as.numeric(sortedGini),
                    stringsAsFactors = FALSE)

# 保留前15个基因进行展示（如果基因数量不足15，保留全部）
nTop <- min(15, nrow(impDF))
impDF <- impDF[1:nTop,]

# 将基因名转换为因子，并按重要性逆序排列
impDF$Gene <- factor(impDF$Gene, levels = rev(impDF$Gene))

# 使用ggplot2和viridis绘制气泡图
p1 <- ggplot(impDF, aes(x = Importance, y = Gene)) +
  geom_segment(aes(x = 0, xend = Importance, y = Gene, yend = Gene), color = "grey80", size = 1.2) +
  geom_point(aes(color = Importance, size = Importance), shape = 21, fill = "white", stroke = 1.2) +
  scale_color_viridis(option = "plasma", direction = -1) +
  scale_size_continuous(range = c(3,7)) +
  labs(x = "Importance", y = "", title = "Gene Importance") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# 保存气泡图到PDF文件
pdf(file = "Gene_Importance_Bubble.pdf", width = 7, height = 5)
print(p1)
dev.off()
# 构建所有基因评分的数据框
geneScores <- data.frame(Gene = names(sortedGini), Importance = as.numeric(sortedGini),
                         stringsAsFactors = FALSE)

# 将所有基因评分输出到文本文件（例如CSV或制表符分隔文件）
write.table(geneScores, file = "All_Gene_Scores.txt", sep = "\t", quote = FALSE, row.names = FALSE)
