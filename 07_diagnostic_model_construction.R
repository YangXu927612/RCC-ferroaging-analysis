# ==============================================================================
# 07_diagnostic_model_construction.R
# One-stop integrated R script for gene analysis pipeline.
# Builds a multi-gene diagnostic model (logistic regression / nomogram),
# generates calibration curves, decision curves, and ROC/AUC plots.
# ==============================================================================

# ============================ Load required packages =============================
library(rms)
library(regplot)
library(rmda)
library(pROC)
library(tools)

# Step 0: Set the working directory (relative path; no absolute paths used).
working_directory <- getwd()
if (!dir.exists(working_directory)) {
  stop(paste("Working directory does not exist: ", working_directory))
}
setwd(working_directory)

# Step 1: 进度条与分析步骤
cat("==== 基因分析pipeline启动 ====\n")
steps <- c(
  "创建输出目录", "读取表达数据", "读取重要基因表", "筛选基因",
  "样本结构化", "特征准备", "建模/列线图", "校准曲线", "决策曲线"
)
pb <- txtProgressBar(min=0, max=length(steps), width=40, style=3)

# Step 2: 创建输出文件夹
result_dir <- "Analysis_Result"
if (!dir.exists(result_dir)) dir.create(result_dir)
setwd(result_dir)
setTxtProgressBar(pb, 1)

# Step 3: 读取数据
input_expr <- file.path("..", "merged_after_batch_removal.csv")
gene_file  <- file.path("..", "Intersect_Genes.csv")
if (!file.exists(input_expr)) stop(sprintf("未找到文件: %s", input_expr))
if (!file.exists(gene_file)) stop(sprintf("未找到文件: %s", gene_file))
dat_expr <- tryCatch(
  read.table(input_expr, sep = ",", header = TRUE, row.names = 1, check.names = FALSE),
  error = function(e) stop("表达文件读取失败")
)
setTxtProgressBar(pb, 2)

# 读取基因列表（CSV格式，有表头Gene）
gene_tbl <- tryCatch(
  read.csv(gene_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) stop("基因列表读取失败")
)
if (ncol(gene_tbl)<1) stop("基因文件无有效列")
feature_genes <- as.character(gene_tbl[,1])
feature_genes <- trimws(feature_genes)
feature_genes <- feature_genes[feature_genes != ""]
if (length(feature_genes)<2) stop("特征基因数过少")
cat(sprintf("读取到 %d 个诊断基因\n", length(feature_genes)))
setTxtProgressBar(pb, 3)

# Step 4: 筛选表达矩阵
dat_expr <- dat_expr[gsub("-", "_", rownames(dat_expr)), , drop=FALSE] # 格式统一
found_genes <- feature_genes %in% rownames(dat_expr)
if (sum(found_genes)==0) stop("IntersectionGenes.csv中基因在表达数据中均找不到")
if (any(!found_genes)) warning(paste0("以下基因不在表达矩阵中：", paste(feature_genes[!found_genes],collapse=",")))
dat_expr_filt <- dat_expr[feature_genes[found_genes],,drop=FALSE]
setTxtProgressBar(pb, 4)

# Step 5: 转置+分组变量赋值
df_expr <- as.data.frame(t(dat_expr_filt))
sample_names <- rownames(df_expr)
groups <- gsub("(.*)_([A-Za-z0-9]+)$", "\\2", sample_names)
df_expr$GroupType <- groups
if (length(unique(groups)) < 2) warning("分组数不足二，分析或有问题")
setTxtProgressBar(pb, 5)

# Step 6: 建模环境准备
ddinfo <- datadist(df_expr)
options(datadist="ddinfo")
setTxtProgressBar(pb, 6)

# Step 7: 组装回归公式
model_vars <- setdiff(colnames(df_expr), "GroupType")
reg_formula <- as.formula(paste("GroupType ~", paste(model_vars, collapse=" + ")))
setTxtProgressBar(pb, 7)

# Step 8: 模型拟合、列线图绘制
lrm_fit <- lrm(reg_formula, data=df_expr, x=TRUE, y=TRUE)

# 提取模型系数并显示
cat("\n==== 模型系数解释 ====\n")
coef_summary <- summary(lrm_fit)
cat("各基因的模型系数方向：\n")
for (i in 1:nrow(coef_summary)) {
  gene_name <- rownames(coef_summary)[i]
  coef_val <- coef_summary[i, "Effect"]
  if (!is.na(coef_val)) {
    direction <- ifelse(coef_val > 0, "正相关(↑表达→↑风险)", "负相关(↓表达→↑风险)")
    cat(sprintf("  %s: %.4f %s\n", gene_name, coef_val, direction))
  }
}
cat("====\n\n")

nomo_obj <- nomogram(
  lrm_fit, fun = plogis,
  fun.at = c(0.001,0.1,0.3,0.5,0.7,0.9,0.99), lp=FALSE, funlabel="Disease Risk"
)
pdf("Nomogram_Plot.pdf", width=11, height=6)
plot(nomo_obj)
dev.off()
cat(">> 列线图输出至: Nomogram_Plot.pdf\n")
cat(">> 注：列线图中基因的排序方向反映模型系数方向，不是表达水平本身\n")
setTxtProgressBar(pb, 8)

# Step 9: 校准曲线
calibrate_obj <- calibrate(lrm_fit, method="boot", B=800)
pdf("Calibration_Curve.pdf", width=5.5, height=5.5)

# 手动绘制校准曲线（不含C.L.）
cal_data <- calibrate_obj[, c("predy", "calibrated.orig", "calibrated.corrected")]
plot(cal_data[,"predy"], cal_data[,"calibrated.corrected"],
     type="l", lwd=2, col="black",
     xlab="Predicted", ylab="Observed",
     xlim=c(0,1), ylim=c(0,1))
# 添加理想线（对角线）
abline(a=0, b=1, lty=2, col="gray50")
# 添加原始校准线
lines(cal_data[,"predy"], cal_data[,"calibrated.orig"], lty=1, lwd=1, col="darkblue")
# 添加图例
legend("bottomright", legend=c("Ideal", "Apparent", "Bias-corrected"),
       lty=c(2,1,1), lwd=c(1,1,2), col=c("gray50","darkblue","black"), bty="n")

dev.off()
cat(">> 校准曲线输出至: Calibration_Curve.pdf\n")

# ==== Step 10: 决策曲线分析 (DCA) - 联合模型 + 单基因 ====
df_expr$GroupType <- as.numeric(factor(df_expr$GroupType)) - 1
if(!all(df_expr$GroupType %in% c(0, 1))) stop("GroupType不是严格的0/1!")
set.seed(123)

# 构建DCA数据框
df_expr_dca <- df_expr
df_expr_dca$Diagnostic_Model <- predict(lrm_fit, newdata=df_expr, type="fitted")

# 联合模型DCA
dca_combined <- decision_curve(
  GroupType ~ Diagnostic_Model,
  data = df_expr_dca,
  thresholds = seq(0, 1, by = 0.01),
  family = binomial(link = "logit"),
  bootstraps = 100
)

# 每个单基因分别建DCA对象
dca_gene_list <- list()
for (g in model_vars) {
  gene_formula <- as.formula(paste("GroupType ~", g))
  dca_gene_list[[g]] <- tryCatch(
    decision_curve(
      gene_formula,
      data = df_expr_dca,
      thresholds = seq(0, 1, by = 0.01),
      family = binomial(link = "logit"),
      bootstraps = 100
    ),
    error = function(e) {
      cat(sprintf("  警告: 基因 %s 的DCA计算失败: %s\n", g, e$message))
      NULL
    }
  )
}
# 去除失败的
dca_gene_list <- dca_gene_list[!sapply(dca_gene_list, is.null)]

# 合并所有DCA对象为列表
all_dca <- c(list(dca_combined), dca_gene_list)

# DCA配色：联合模型红色，单基因依次不同颜色
dca_colors <- c("red", "deepskyblue", "forestgreen", "orange",
                "purple", "magenta", "gold", "brown", "gray40", "cyan")
while(length(dca_colors) < length(all_dca)) {
  dca_colors <- c(dca_colors, rainbow(length(all_dca) - length(dca_colors)))
}
dca_colors <- dca_colors[1:length(all_dca)]

# 图例标签
dca_labels <- c("Diagnostic_Model", names(dca_gene_list))

pdf(file = "DCA.pdf", width = 9, height = 8)
# 用layout分割：上方图例区(1)，下方主图区(2)
layout(matrix(c(1, 2), nrow = 2), heights = c(0.8, 5))

# --- 上方：单独画图例 ---
par(mar = c(0, 4, 0.5, 2))
plot.new()
legend_labels <- c(dca_labels, "Treat All", "Treat None")
legend_colors <- c(dca_colors, "black", "black")
legend_lty <- c(rep(1, length(dca_labels)), 1, 2)
legend("center",
       legend = legend_labels,
       col = legend_colors,
       lty = legend_lty,
       lwd = 2,
       ncol = 3,
       bty = "n",
       cex = 1.15)

# --- 下方：主图 ---
par(mar = c(5, 4, 0.5, 2))
plot_decision_curve(
  all_dca,
  xlab = "Threshold Probability",
  col = dca_colors,
  confidence.intervals = FALSE,
  standardize = TRUE,
  cost.benefit.axis = TRUE,
  legend.position = "none"
)
dev.off()

cat(">> DCA决策曲线（联合模型+单基因）输出至 DCA.pdf\n")




# Step 11: 输出主要分析表格
write.csv(df_expr, file="Filtered_Expression_Matrix.csv", quote=F)
write.csv(as.data.frame(dca_combined$derived.data), file="Decision_Curve_Data.csv")
pred_probs <- predict(lrm_fit, newdata=df_expr, type="fitted")
df_out <- data.frame(Sample=rownames(df_expr), GroupType=df_expr$GroupType, Predicted_Prob=pred_probs)
write.csv(df_out, file="Sample_Predicted_Probabilities.csv")

# Step 12: 模型系数/OR/CI结果表
coef_df <- as.data.frame(summary(lrm_fit))
if("Effect" %in% names(coef_df) & "S.E." %in% names(coef_df)) {
  coef_df$Effect <- as.numeric(as.character(coef_df$Effect))
  coef_df$`S.E.` <- as.numeric(as.character(coef_df$`S.E.`))
  coef_df$OR <- exp(coef_df$Effect)
  coef_df$OR_low <- exp(coef_df$Effect - 1.96 * coef_df$`S.E.`)
  coef_df$OR_high <- exp(coef_df$Effect + 1.96 * coef_df$`S.E.`)
  write.csv(coef_df, file="Model_Coefficients.csv", row.names=FALSE)
  if ("Lower 0.95" %in% names(coef_df) & "Upper 0.95" %in% names(coef_df)) {
    write.csv(coef_df, file="Model_Coefficients_With_CI.csv", row.names=FALSE)
  }
  or_selected <- coef_df[!is.na(coef_df$OR) & coef_df$OR > 1, ]
  write.csv(or_selected, file="Model_OR_GT1.csv", row.names=FALSE)
  if("P" %in% names(coef_df)) {
    sig_coef <- coef_df[!is.na(coef_df$P) & coef_df$P < 0.05, ]
    write.csv(sig_coef, file="Model_Significant_Coefficients.csv", row.names=FALSE)
  }
}

# Step 13：标注列线图（regplot风格，带分布条、P值、观测标注）
# 注意：regplot 无法用 pdf()/png() 包裹导出，需在 RStudio Plot 面板中手动 Export 保存
median_point <- sapply(df_expr[,model_vars,drop=FALSE], median)
median_df <- as.data.frame(t(median_point)); median_df$GroupType <- 1
rownames(median_df) <- "AllMedian"

cat(">> 正在绘制 regplot 列线图，请在 RStudio 右下角 Plots 面板中手动 Export 保存为 PDF/PNG\n")
regplot(
  lrm_fit, showP=TRUE, rank="sd", distribution=TRUE,
  observation=median_df, title="Prediction Nomogram"
)
cat(">> regplot 列线图已在 Plots 面板展示，请手动保存\n")
# ==== 联合模型 ROC（主模型红色，风格与单基因统一）====


roc_y <- as.numeric(df_expr$GroupType)
roc_pred <- pred_probs
roc_obj <- roc(roc_y, roc_pred, levels = c(0,1), direction = "<")

# 主模型/单基因统一色板（前几个色都与后面单基因一致）
my_cols <- c("red", "deepskyblue", "forestgreen", "orange", 
             "purple", "gray40", "black", "magenta", "gold", "brown")
while(length(my_cols) < length(model_vars) + 1) { 
  my_cols <- c(my_cols, rainbow(length(model_vars) + 1 - length(my_cols)))
}

## ====== 联合主模型ROC ======
pdf("CombinedGenes_ROC_SCI.pdf", width = 6, height = 6)
par(mar = c(5,6,4,2)+0.1, cex = 1.3)
plot(1-roc_obj$specificities, roc_obj$sensitivities, type = "l",
     col = my_cols[1], lwd = 4, lty = 1,
     xlab = expression("1 - Specificity"), ylab = "Sensitivity",
     main = " Model ROC Curve", 
     cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45,
     xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2, col = "gray70", lwd = 2)
auc_val <- as.numeric(auc(roc_obj))
legend("bottomright", legend = sprintf("AUC = %.3f", auc_val), col = my_cols[1], 
       lwd = 4, lty=1, bty = "n", cex = 1.2)
dev.off()
cat(">> 合并主模型 ROC 曲线 (红色) 已输出至 CombinedGenes_ROC_SCI.pdf\n")
# ==== 单基因多色ROC优化版 ====
roc_list <- list(); auc_list <- c(); leglab <- c(); auc_dir <- c()
group_names <- levels(factor(roc_y))

pdf("IndividualGenes_ROC_SCI.pdf", width = 8, height = 8)
par(mar = c(5, 6, 4, 2)+0.1, cex = 1.3)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     xlab = expression("1 - Specificity"), ylab = "Sensitivity",
     main = "Gene ROC Curves", 
     cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45
)
abline(0, 1, lty = 2, col = "gray70", lwd = 2)

for(i in seq_along(model_vars)) {
  g <- model_vars[i]
  # 先用direction="auto"自动判断
  cur_roc <- roc(roc_y, df_expr[[g]], levels = c(0,1), direction = "auto")
  cur_auc <- as.numeric(auc(cur_roc))
  
  # 为确保AUC反映“最强判别能力”，若AUC<0.5, 则取1-AUC, 并标记方向
  flip <- FALSE
  if(cur_auc < 0.5) {
    flip <- TRUE
    cur_auc <- 1 - cur_auc  # 修正为正向（高表达区分度而不是方向差，见下）
    # 必须重新画ROC曲线，用反向direction
    cur_roc <- roc(roc_y, df_expr[[g]], levels = c(0,1), direction = ifelse(cur_roc$direction==">", "<", ">"))
  }
  
  lines(1-cur_roc$specificities, cur_roc$sensitivities, col = my_cols[i+1], lwd = 3)
  roc_list[[g]] <- cur_roc
  auc_list <- c(auc_list, cur_auc)
  
  # 写入方向注释
  mean0 <- mean(df_expr[[g]][roc_y==0], na.rm=TRUE)
  mean1 <- mean(df_expr[[g]][roc_y==1], na.rm=TRUE)
  direction_note <- ifelse(mean1 > mean0, paste0("(组1高)"), paste0("(组0高)"))
  auc_dir <- c(auc_dir, direction_note)
  
  # 图例写明方向和AUC
  leglab <- c(leglab, sprintf("%s  AUC=%.3f", g, cur_auc))
  
}
legend("bottomright", legend = leglab, col = my_cols[2:(length(model_vars)+1)], 
       lwd = 3, lty = 1, bty = "n", cex = 1.15, y.intersp = 1.18)
dev.off()

# AUC同步输出加强版
gene_auc_df <- data.frame(
  Gene = model_vars,
  AUC = auc_list,
  Direction = auc_dir
)
write.csv(gene_auc_df, file = "IndividualGenes_AUC.csv", row.names = FALSE)
cat(">> 单基因 ROC 曲线和 AUC 已输出至 IndividualGenes_ROC_SCI.pdf / IndividualGenes_AUC.csv\n")

# ==== 新增：诊断模型基因的差异箱线图（Nature风格）====
cat(">> 生成诊断模型基因的差异箱线图...\n")

# 准备数据用于箱线图
df_boxplot <- data.frame()
for (gene in model_vars) {
  temp_df <- data.frame(
    Gene = gene,
    Expression = df_expr[[gene]],
    Group = factor(df_expr$GroupType, labels = c("Control", "Disease"))
  )
  df_boxplot <- rbind(df_boxplot, temp_df)
}

# 计算p值（t检验）
p_values <- c()
for (gene in model_vars) {
  group0 <- df_expr[[gene]][df_expr$GroupType == 0]
  group1 <- df_expr[[gene]][df_expr$GroupType == 1]
  p_val <- t.test(group0, group1)$p.value
  p_values <- c(p_values, p_val)
}

# 创建p值标签
p_labels <- sapply(p_values, function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
})

# 绘制Nature风格的箱线图
pdf("Diagnostic_Genes_Boxplot.pdf", width = 6, height = 5)

# 使用ggplot2绘制
library(ggplot2)
library(ggpubr)

# 计算y轴最大值用于添加显著性标记
y_max <- max(df_boxplot$Expression, na.rm = TRUE)
y_min <- min(df_boxplot$Expression, na.rm = TRUE)
y_range <- y_max - y_min

# 创建箱线图
p <- ggplot(df_boxplot, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = 21,
    outlier.size = 2,
    outlier.fill = "white",
    outlier.stroke = 1,
    alpha = 0.8,
    lwd = 0.8
  ) +
  scale_fill_manual(
    values = c("Control" = "#E8E8E8", "Disease" = "#FF6B6B"),
    name = "Group"
  ) +
  scale_color_manual(
    values = c("Control" = "#333333", "Disease" = "#CC0000")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
  ) +
  labs(
    y = "Expression Level",
    title = "Diagnostic Model Genes Expression Comparison"
  )

# 添加显著性标记
y_pos <- y_max + y_range * 0.05
for (i in seq_along(model_vars)) {
  p <- p + annotate(
    "text",
    x = i,
    y = y_pos,
    label = p_labels[i],
    size = 5,
    fontface = "bold"
  )
}

print(p)
dev.off()

cat(">> 差异箱线图已输出至 Diagnostic_Genes_Boxplot.pdf\n")

# 输出p值表格
p_value_df <- data.frame(
  Gene = model_vars,
  P_Value = p_values,
  Significance = p_labels
)
write.csv(p_value_df, file = "Diagnostic_Genes_PValues.csv", row.names = FALSE)
cat(">> p值表格已输出至 Diagnostic_Genes_PValues.csv\n")


# ==== 新增：多因素Logistic回归 OR值森林图 ====
cat(">> 计算OR值并绘制森林图...\n")

# 用glm拟合多因素logistic回归（与lrm模型一致）
glm_fit <- glm(reg_formula, data = df_expr, family = binomial(link = "logit"))
glm_summary <- summary(glm_fit)

# 提取系数（排除截距）
coef_table <- coef(glm_summary)
coef_table <- coef_table[rownames(coef_table) != "(Intercept)", , drop = FALSE]

# 计算OR及95%CI
or_df <- data.frame(
  Gene     = rownames(coef_table),
  Beta     = coef_table[, "Estimate"],
  SE       = coef_table[, "Std. Error"],
  Z        = coef_table[, "z value"],
  P        = coef_table[, "Pr(>|z|)"]
)
or_df$OR       <- exp(or_df$Beta)
or_df$OR_lower <- exp(or_df$Beta - 1.96 * or_df$SE)
or_df$OR_upper <- exp(or_df$Beta + 1.96 * or_df$SE)

# P值星号标记
or_df$Sig <- ifelse(or_df$P < 0.001, "***",
             ifelse(or_df$P < 0.01, "**",
             ifelse(or_df$P < 0.05, "*", "")))

# OR标签文本
or_df$Label <- sprintf("%.2f (%.2f-%.2f)%s", or_df$OR, or_df$OR_lower, or_df$OR_upper, or_df$Sig)

# 按OR值排序（从大到小）
or_df <- or_df[order(or_df$OR), ]
or_df$Gene <- factor(or_df$Gene, levels = or_df$Gene)

# 输出OR表格
write.csv(or_df, file = "Logistic_OR_ForestPlot_Data.csv", row.names = FALSE)
cat(">> OR值表格已输出至 Logistic_OR_ForestPlot_Data.csv\n")

# 绘制森林图
library(ggplot2)

# 颜色：OR>1红色（风险），OR<1蓝色（保护）
or_df$Direction <- ifelse(or_df$OR > 1, "Risk", "Protective")

pdf("OR_ForestPlot.pdf", width = 8, height = max(4, nrow(or_df) * 0.8 + 1.5))
p_forest <- ggplot(or_df, aes(x = OR, y = Gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", size = 0.8) +
  geom_errorbar(aes(xmin = OR_lower, xmax = OR_upper),
                width = 0.25, linewidth = 0.8, color = "grey30", orientation = "y") +
  geom_point(aes(fill = Direction), shape = 21, size = 4.5, stroke = 0.6, color = "black") +
  scale_fill_manual(values = c("Risk" = "#BC3C29", "Protective" = "#0072B5"),
                    name = "Direction") +
  scale_x_log10() +
  geom_text(aes(x = OR, label = Label),
            vjust = 2.2, hjust = 0.5, size = 3.8, color = "black") +
  theme_bw(base_size = 14) +
  labs(x = "Odds Ratio (95% CI, log scale)", y = NULL,
       title = " ") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(face = "bold", size = 13),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.margin = ggplot2::margin(10, 10, 10, 10, unit = "pt")
  )
print(p_forest)
dev.off()
cat(">> OR森林图已输出至 OR_ForestPlot.pdf\n")

# 收尾
cat("==== 全流程执行完毕！====\n")
