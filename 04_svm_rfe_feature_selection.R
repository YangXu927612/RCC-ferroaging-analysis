
  # ==============================================================================
# 04_svm_rfe_feature_selection.R
# SVM-RFE (Support Vector Machine - Recursive Feature Elimination) for
# candidate gene selection from the ferro-aging-related differentially
# expressed genes.
# ==============================================================================

library(e1071)                                     # load e1071 if installed


# Set random seed for reproducibility.
set.seed(12345)                                      # fixed seed

# Define input file paths (relative paths only; no absolute paths used).
workDir <- getwd()                                       # use current directory
inputFile <- file.path(workDir, "merged_after_batch_removal.csv")  # expression matrix

# 自动识别 Intersection_ 前缀的txt文件作为基因列表
gene_list_files <- list.files(workDir, pattern = "^Intersection_.*\\.txt$", full.names = TRUE)
if (length(gene_list_files) == 0) {
  stop("未在工作目录下找到 Intersection_ 前缀的txt文件！")
}
geneListFile <- gene_list_files[1]  # 使用第一个匹配的文件
message(sprintf("自动识别基因列表文件: %s", basename(geneListFile)))

# 设置当前工作目录
setwd(workDir)                                       # 设置工作目录

# 判断输入文件是否存在，不存在则停止执行
if (!file.exists(inputFile)) {                      # 检查表达数据文件是否存在
  stop("表达数据文件不存在，请检查文件路径！")         # 若不存在则报错停止执行
}

if (!file.exists(geneListFile)) {                   # 检查基因列表文件是否存在
  stop("基因列表文件不存在，请检查文件路径！")         # 若不存在则报错停止执行
}

# 读取要分析的基因列表（自动判断是否有表头）
message("读取基因列表文件...")
geneListRaw <- readLines(geneListFile)
geneListRaw <- trimws(geneListRaw)
geneListRaw <- geneListRaw[geneListRaw != ""]
# 如果第一行是 "Gene" 或 "Genes" 等表头，则跳过
if (tolower(geneListRaw[1]) %in% c("gene", "genes", "genename", "gene_name")) {
  geneList <- geneListRaw[-1]
} else {
  geneList <- geneListRaw
}
geneList <- unique(geneList)
message(sprintf("共读取 %d 个基因", length(geneList)))

# ---------------------- 定义辅助函数 ----------------------

# 定义获取SVM权重的函数（用于计算特征重要性）
getWeights <- function(test.fold, X) {
  # 根据是否指定测试折，选择训练数据
  train.data <- X
  if (!is.null(test.fold)) {                         # 如果提供测试折编号
    train.data <- X[-test.fold, ]                    # 排除测试折数据作为训练数据
  }
  # 使用线性核训练SVM模型，cost设为10，关闭内部数据缩放
  svmModel <- svm(train.data[, -1], train.data[, 1],
                  cost = 10, cachesize = 500,
                  scale = FALSE, type = "C-classification",
                  kernel = "linear")                # 训练线性SVM分类器
  # 计算并返回权重向量（模型系数乘支持向量）
  return(t(svmModel$coefs) %*% svmModel$SV)          # 返回SVM权重向量
}

# 定义SVM-RFE算法函数，用于特征选择
svmRFE <- function(X, k = 1, halve.above = 5000) {
  totalFeatures <- ncol(X) - 1                       # 计算除分组列外的特征总数
  cat("总特征数量：", totalFeatures, "\n")           # 输出总特征数量
  cat("开始数据标准化...\n")                         # 输出步骤提示
  
  # 对数据进行标准化处理（不对分组信息处理）
  X[, -1] <- scale(X[, -1])                          # 标准化特征数据
  cat("数据标准化完成！\n")                           # 输出标准化完成信息
  
  # 初始化进度条，显示特征递归过程
  pb <- txtProgressBar(min = 1, max = totalFeatures, style = 3)  # 初始化进度条
  
  survivingIndices <- 1:totalFeatures                # 初始化存活特征索引
  rankingPosition  <- totalFeatures                  # 初始化排名位置
  rankedFeatures <- vector(length = totalFeatures)   # 初始化最终特征排名列表
  
  # 循环递归消除特征，直至所有特征均被排序
  while (length(survivingIndices) > 0) {
    # 判断是否进行多重SVM-RFE
    if (k > 1) {                                     # 如果k大于1则使用多重验证方式
      # 生成随机折，保证每个折均匀分布
      folds <- rep(1:k, length.out = nrow(X))
      folds <- sample(folds)                         # 打乱折序
      foldsList <- lapply(1:k, function(x) which(folds == x))  # 将折编号分组
      
      # 对每个折进行训练，获取权重向量
      weightList <- lapply(foldsList, getWeights, X[, c(1, 1 + survivingIndices)])
      weightMatrix <- do.call(rbind, weightList)      # 合并各折的权重向量
      
      # 标准化每个权重向量
      weightMatrix <- t(apply(weightMatrix, 1, function(w) w / sqrt(sum(w^2) + 1e-6)))
      
      # 计算每个特征的评分，利用均值和标准差
      v <- weightMatrix^2                           # 计算平方值
      vBar <- apply(v, 2, mean)                       # 计算均值
      vSD <- apply(v, 2, sd) + 1e-6                   # 计算标准差，防止除零
      scores <- vBar / vSD                            # 得到每个特征的评分
    } else {                                         # 若只进行单次SVM-RFE
      weights <- getWeights(NULL, X[, c(1, 1 + survivingIndices)])  # 计算权重向量
      scores <- weights^2                           # 得到评分（平方）
    }
    
    # 对当前存活特征按评分进行排序（从小到大）
    ranking <- sort(scores, index.return = TRUE)$ix   # 获取排序索引
    if (length(survivingIndices) == 1) {             # 特征仅剩一时
      ranking <- 1                                 # 设置索引为1
    }
    
    # 判断是否需要一次性剔除多个特征（特征数大于阈值时）
    if (length(survivingIndices) > halve.above) {   # 如果存活特征数量大于阈值
      numFeatures <- length(survivingIndices)        # 当前存活特征数
      numToCut <- round(numFeatures / 2)             # 一次剔除一半特征
      cat("当前特征数:", numFeatures, "，剔除", numToCut, "个特征...\n")  # 提示剔除信息
      # 重设进度条范围
      pb <- txtProgressBar(min = 1, max = numFeatures - numToCut, style = 3)
    } else {
      numToCut <- 1                                # 否则每次剔除一个特征
    }
    
    # 将剔除的特征保存至排名列表中
    rankedFeatures[rankingPosition:(rankingPosition - numToCut + 1)] <- 
      survivingIndices[ranking[1:numToCut]]
    rankingPosition <- rankingPosition - numToCut   # 更新排名位置
    survivingIndices <- survivingIndices[-ranking[1:numToCut]]  # 移除已剔除特征
    
    # 更新并显示进度条
    setTxtProgressBar(pb, length(survivingIndices))  # 更新进度条显示
    flush.console()                                # 刷新控制台
  }
  
  # 关闭进度条
  close(pb)                                          # 结束进度条显示
  
  # 返回特征的排名结果
  return(rankedFeatures)                             # 返回最终特征排名向量
}

# 包装SVM-RFE算法的函数，用于交叉验证折内的特征选择
svmRFE_wrap <- function(test.fold, X, ...) {
  # 分离训练集与测试集
  trainingData <- X[-test.fold, ]                   # 训练数据
  testingData  <- X[test.fold, ]                    # 测试数据
  
  # 调用svmRFE函数对训练集进行特征排序
  featuresRanked <- svmRFE(trainingData, ...)
  
  # 返回包含特征排名以及训练集、测试集的行号信息的列表
  return(list(feature_ids = featuresRanked,
              train_ids = row.names(trainingData),
              test_ids = row.names(testingData)))
}

# 定义特征排序结果写入函数
WriteFeatures <- function(results, input, save = TRUE, file = 'features_ranked.txt') {
  # 对各折的特征排序结果求平均排名
  featureID <- sort(apply(sapply(results, function(x) {
    # 此处增加冗余判断，防止出现空数据
    if(is.null(x$feature)) return(rep(NA, ncol(input) - 1))
    sort(x$feature, index.return = TRUE)$ix
  }), 1, mean, na.rm = TRUE), index.return = TRUE)$ix
  avgRank <- sort(apply(sapply(results, function(x) {
    if(is.null(x$feature)) return(rep(NA, ncol(input) - 1))
    sort(x$feature, index.return = TRUE)$ix
  }), 1, mean, na.rm = TRUE), index.return = TRUE)$x
  
  # 获取特征名称（排除第一列分组信息）
  featureNames <- colnames(input[, -1])[featureID]
  # 构建数据框存储排序结果
  featuresRanked <- data.frame(FeatureName = featureNames,
                               FeatureID = featureID,
                               AvgRank = avgRank)
  # 根据save参数决定是否写入文件
  if (save) {
    write.table(featuresRanked, file = file, quote = FALSE, row.names = FALSE)
  } else {
    return(featuresRanked)
  }
}

# 定义泛化误差验证函数，对给定特征数量进行评估
FeatSweep_wrap <- function(i, results, input) {
  # 对每个交叉验证折进行SVM调优，并获取模型性能
  svmList <- lapply(results, function(x) {
    # 使用tune函数寻找最优超参数，包含gamma与cost的搜索
    tuned <- tune(svm,
                  train.x = input[x$train_ids, 1 + x$feature_ids[1:i]],
                  train.y = input[x$train_ids, 1],
                  ranges = list(gamma = 2^(-12:0), cost = 2^(-6:6)),
                  tunecontrol = tune.control(sampling = "fix"))
    # 重新用最优参数进行模型验证
    model_perf <- tune(svm,
                       train.x = input[x$train_ids, 1 + x$feature_ids[1:i]],
                       train.y = input[x$train_ids, 1],
                       validation.x = input[x$test_ids, 1 + x$feature_ids[1:i]],
                       validation.y = input[x$test_ids, 1],
                       ranges = tuned$best.parameters,
                       tunecontrol = tune.control(sampling = "fix"))$perf
    return(model_perf)
  })
  # 计算所有折的平均错误率
  errorMean <- mean(sapply(svmList, function(perf) perf$error))
  return(list(svm_list = svmList, error = errorMean))
}

# 包装泛化误差验证函数，增加进度条和剩余时间估算
FeatSweep_wrap_progress <- function(i, results, input) {
  startTime <- Sys.time()                           # 记录开始时间
  # 执行泛化误差评估
  result <- FeatSweep_wrap(i, results, input)
  elapsedTime <- as.numeric(Sys.time() - startTime, units = "secs")  # 计算已用时间
  
  # 估算剩余时间（假设评估总共30个特征）
  remainingTime <- if (i > 0) elapsedTime * (30 - i) / i else NA
  # 输出当前进度信息
  cat(sprintf("已完成特征数：%d/30, 耗时：%.2f秒, 预计剩余时间：%.2f秒\n",
              i, elapsedTime, remainingTime))
  return(result)
}

# 定义绘图函数，绘制泛化误差曲线
# 参数 optimal_n: 可选，指定最优特征数量，确保与准确率图一致
PlotErrors <- function(errors, errors2 = NULL, no.info = 0.5,
                       ylim = range(c(errors, errors2), na.rm = TRUE),
                       xlab = 'Number of Features', ylab = '10 x CV Error',
                       optimal_n = NULL) {
  oldPar <- par(mar = c(5, 5, 4, 5), xpd = TRUE)    # 保存原始图形参数

  # 内部函数，添加数据曲线与标记
  AddLine <- function(x, col = '#4B4453', opt_n = NULL) {
    lines(which(!is.na(x)), na.omit(x), col = col, lwd = 3)  # 绘制曲线
    # 如果提供了最优特征数，使用它；否则自动计算
    x_min <- if (!is.null(opt_n)) opt_n else which.min(x)
    y_min <- x[x_min]                               # 获取对应的错误率
    points(x_min, y_min, col = '#C34A36', pch = 19, cex = 1.5)  # 标记最小值
    text(x_min, y_min, labels = paste0("n=", x_min, "\n(", sprintf("%.4f", y_min), ")"),
         pos = 2, col = 'darkred', cex = 1.2, font = 2)  # 添加文字说明
  }

  # 绘制空图框
  plot(errors, type = 'n', ylim = ylim, xlab = xlab, ylab = ylab,
       main = "Generalization Error vs Number of Features",
       col.main = "darkblue", col.lab = "darkblue", cex.lab = 1.2, cex.axis = 1.1)
  grid(col = "lightgray", lty = "dotted")          # 添加网格线
  AddLine(errors, opt_n = optimal_n)               # 绘制错误率曲线
  if (!is.null(errors2)) {                         # 如果有第二组错误率数据
    AddLine(errors2, col = 'forestgreen')           # 绘制第二条曲线
  }
  abline(h = no.info, lty = 2, col = 'darkgray')     # 添加无信息判别线
  on.exit(par(oldPar))                              # 恢复原始图形参数
}

# 定义绘图函数，绘制准确率曲线（准确率 = 1 - 错误率）
# 参数 optimal_n: 可选，指定最优特征数量，确保与错误率图一致
PlotAccuracy <- function(accuracy, accuracy2 = NULL, no.info = 0.5,
                         ylim = range(c(accuracy, accuracy2), na.rm = TRUE),
                         xlab = 'Number of Features', ylab = '10 x CV Accuracy',
                         optimal_n = NULL) {
  oldPar <- par(mar = c(5, 5, 4, 5), xpd = TRUE)    # 保存原始图形参数

  # 内部函数，添加数据曲线与标记
  AddLine <- function(x, col = '#4B4453', opt_n = NULL) {
    lines(which(!is.na(x)), na.omit(x), col = col, lwd = 3)  # 绘制曲线
    # 如果提供了最优特征数，使用它；否则自动计算
    x_max <- if (!is.null(opt_n)) opt_n else which.max(x)
    y_max <- x[x_max]                               # 获取对应的准确率值
    points(x_max, y_max, col = 'gold', pch = 19, cex = 1.5)  # 标记最高值
    text(x_max, y_max, labels = paste0("n=", x_max, "\n(", sprintf("%.4f", y_max), ")"),
         pos = 4, offset = 1, col = 'darkred', cex = 1.2, font = 2)  # 添加文字说明
  }

  # 绘制空图框
  plot(accuracy, type = 'n', ylim = ylim, xlab = xlab, ylab = ylab,
       main = "CV Accuracy vs Number of Features",
       col.main = "darkblue", col.lab = "darkblue", cex.lab = 1.2, cex.axis = 1.1)
  grid(col = "lightgray", lty = "dotted")          # 添加网格线
  AddLine(accuracy, opt_n = optimal_n)             # 绘制准确率曲线
  if (!is.null(accuracy2)) {                       # 如果有第二组数据
    AddLine(accuracy2, col = 'forestgreen')         # 绘制第二条曲线
  }
  abline(h = no.info, lty = 2, col = 'darkgray')     # 添加无信息判别线
  on.exit(par(oldPar))                              # 恢复原始图形参数
}

# ---------------------- 数据读取与预处理 ----------------------

# 读取输入数据文件，注意检查分隔符和文件格式
message("读取基因表达数据...")
data <- read.csv(inputFile, header = TRUE, check.names = FALSE, row.names = 1)  # 读取CSV文件

# 判断数据是否为空，若为空则报错
if (nrow(data) == 0) {
  stop("数据为空，请检查输入文件内容！")
}

message(sprintf("原始数据：%d 个基因，%d 个样本", nrow(data), ncol(data)))

# 筛选指定的基因
message("筛选指定的基因...")
# 获取数据中存在的基因（行名中存在的基因）
availableGenes <- intersect(geneList, rownames(data))
message(sprintf("在表达数据中找到 %d 个基因", length(availableGenes)))

if (length(availableGenes) == 0) {
  stop("指定的基因列表中没有在表达数据中找到的基因！请检查基因名称是否匹配。")
}

# 只保留指定的基因
data <- data[availableGenes, , drop = FALSE]
message(sprintf("筛选后数据：%d 个基因，%d 个样本", nrow(data), ncol(data)))

# 转置数据，使样本为行，基因为列
data <- as.data.frame(t(data))                      # 转置数据

# 提取样本分组信息（从行名中提取，格式为"样本_组别"）
message("提取样本分组信息...")
group <- rownames(data)
# 尝试从行名中提取组别信息（假设格式为"样本_组别"）
group <- gsub("(.*)\\_(.*)", "\\2", group)

# 检查是否成功提取组别信息
uniqueGroups <- unique(group)
message(sprintf("检测到 %d 个分组：%s", length(uniqueGroups), paste(uniqueGroups, collapse = ", ")))

# 将分组信息添加到数据中，并转换为因子类型
data <- cbind(group, data)                          # 合并分组信息
colnames(data)[1] <- "group"                        # 设置第一列列名

# 尝试自动识别分组标签，如果无法识别则使用原始标签
if (length(uniqueGroups) == 2) {
  # 如果有两个分组，尝试识别对照组和处理组
  if ("con" %in% uniqueGroups && "tre" %in% uniqueGroups) {
    data$group <- factor(data$group, levels = c("con", "tre"))
    message("分组标签：对照组(con) vs 处理组(tre)")
  } else if ("control" %in% tolower(uniqueGroups) || "normal" %in% tolower(uniqueGroups)) {
    # 如果包含control或normal，将其作为第一个水平
    controlLabel <- uniqueGroups[grep("control|normal", tolower(uniqueGroups))][1]
    treatmentLabel <- setdiff(uniqueGroups, controlLabel)
    data$group <- factor(data$group, levels = c(controlLabel, treatmentLabel))
    message(sprintf("分组标签：%s vs %s", controlLabel, treatmentLabel))
  } else {
    # 否则按字母顺序排列
    data$group <- factor(data$group, levels = sort(uniqueGroups))
    message(sprintf("分组标签：%s", paste(sort(uniqueGroups), collapse = " vs ")))
  }
} else {
  data$group <- factor(data$group)
  message(sprintf("警告：检测到 %d 个分组，SVM-RFE 最适合二分类问题", length(uniqueGroups)))
}

# ---------------------- 运行 SVM-RFE 并获取结果 ----------------------

# 执行SVM-RFE算法进行特征选择，参数：多重折数k=10，特征数阈值设为50
message("执行SVM-RFE特征选择...")
svmRFE(data, k = 10, halve.above = 50)              # 调用SVM-RFE函数

# 设置10折交叉验证参数
nfold <- 10                                         # 定义交叉验证折数为10
numSamples <- nrow(data)                            # 样本总数

# 随机生成交叉验证折索引，并分组
folds <- rep(1:nfold, length.out = numSamples)      # 生成折号向量
folds <- sample(folds)                              # 打乱折号顺序
foldsList <- lapply(1:nfold, function(x) which(folds == x))  # 分组存储

# 在每个交叉验证折上运行SVM-RFE，并保存结果到列表中
message("在交叉验证折上运行SVM-RFE...")
results <- lapply(foldsList, svmRFE_wrap, data, k = 10, halve.above = 50)  # 计算每个折的特征排序

# 保存重要特征的排序结果到变量
message("整理特征排序结果...")
topFeatures <- WriteFeatures(results, data, save = FALSE)  # 获取特征排序数据
write.table(topFeatures, file = "feature_svm.txt", sep = "\t", quote = FALSE, row.names = FALSE)  # 写入文件
message("特征排序结果已保存到 feature_svm.txt")

# ---------------------- 泛化误差验证与进度条显示 ----------------------

# 确定要评估的最大特征数（不超过总特征数）
maxFeatures <- min(30, ncol(data) - 1)              # 最多评估30个特征或所有特征
message(sprintf("将评估前 %d 个特征的性能...", maxFeatures))

# 对前N个特征进行泛化误差验证，使用进度条展示进度和剩余时间估计
progressBar <- txtProgressBar(min = 1, max = maxFeatures, style = 3)  # 初始化进度条
featSweep <- lapply(1:maxFeatures, function(i) {
  res <- FeatSweep_wrap_progress(i, results, data)  # 调用泛化误差验证函数
  setTxtProgressBar(progressBar, i)                 # 更新进度条
  return(res)                                       # 返回当前结果
})
close(progressBar)                                  # 关闭进度条

# 计算数据中最小的分布比例，作为无信息判别线
no.info <- min(prop.table(table(data[, 1])))         # 计算无信息错误率

# 计算每个特征数量对应的平均泛化错误率
errors <- sapply(featSweep, function(x) ifelse(is.null(x), NA, x$error))  # 计算平均错误率

# ---------------------- 最优特征选择 ----------------------

# 根据错误率最小的特征数量选择最优特征（在绘图前计算，确保两图一致）
message("确定最优特征...")
optimalFeatureCount <- which.min(errors)            # 获取错误率最小时对应的特征数
optimalError <- min(errors, na.rm = TRUE)          # 获取最小错误率
optimalAccuracy <- 1 - optimalError                # 计算最优准确率

# ---------------------- 绘图输出 ----------------------

# 绘制泛化错误率曲线图，并保存为PDF文件
message("生成泛化错误率曲线图...")
pdf(file = "errors.pdf", width = 5, height = 5)      # 打开PDF设备
PlotErrors(errors, no.info = no.info, optimal_n = optimalFeatureCount)  # 调用绘图函数绘制错误率图
dev.off()                                           # 关闭PDF设备
message("错误率图已保存到 errors.pdf")

# 绘制交叉验证准确率曲线图（准确率 = 1 - 错误率），保存为PDF文件
message("生成交叉验证准确率曲线图...")
pdf(file = "accuracy.pdf", width = 5, height = 5)      # 打开PDF设备
PlotAccuracy(1 - errors, no.info = no.info, optimal_n = optimalFeatureCount)  # 调用绘图函数绘制准确率图
dev.off()                                           # 关闭PDF设备
message("准确率图已保存到 accuracy.pdf")

# 将错误率和准确率图合并显示，并保存为PDF文件
message("生成合并图表...")
pdf(file = "combined_plots.pdf", width = 10, height = 5)  # 打开PDF设备，设置多图布局
par(mfrow = c(1, 2))                                # 设置两图并排显示
PlotErrors(errors, no.info = no.info, optimal_n = optimalFeatureCount)  # 绘制错误率图，使用统一的最优特征数
PlotAccuracy(1 - errors, no.info = no.info, optimal_n = optimalFeatureCount)  # 绘制准确率图，使用统一的最优特征数
dev.off()                                           # 关闭PDF设备
message("合并图表已保存到 combined_plots.pdf")

message(sprintf("最优特征数量：%d", optimalFeatureCount))
message(sprintf("最小错误率：%.4f", optimalError))
message(sprintf("最高准确率：%.4f", optimalAccuracy))

featureGenes <- topFeatures[1:optimalFeatureCount, 1, drop = FALSE]  # 选择最优特征

# 将最优特征写入文本文件，不含列名和行号
write.table(featureGenes, file = "SVM-RFE.gene.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)  # 写入最优特征文件
message("最优特征基因已保存到 SVM-RFE.gene.txt")

# 输出全部步骤完成提示
cat("========================================\n")
cat("所有步骤执行完毕，结果文件已保存！\n")
cat("========================================\n")
cat(sprintf("输入基因数：%d\n", length(geneList)))
cat(sprintf("最终分析基因数：%d\n", ncol(data) - 1))
cat(sprintf("样本总数：%d\n", nrow(data)))
cat(sprintf("最优特征数：%d\n", optimalFeatureCount))
cat(sprintf("最小错误率：%.4f\n", optimalError))
cat(sprintf("最高准确率：%.4f\n", optimalAccuracy))
cat("========================================\n")
cat("输出文件：\n")
cat("  - feature_svm.txt：所有特征的排名结果\n")
cat("  - SVM-RFE.gene.txt：最优特征基因列表\n")
cat("  - errors.pdf：泛化错误率曲线\n")
cat("  - accuracy.pdf：交叉验证准确率曲线\n")
cat("  - combined_plots.pdf：错误率和准确率合并图\n")
cat("========================================\n")
