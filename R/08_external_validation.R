# ====================================================================
# 08_external_validation.R
# External validation of the multi-gene diagnostic model.
#   - Rebuild the model on the training set (merged_after_batch_removal.csv)
#   - Plot ROC curves on the external validation cohort GSE53757
#   - Also draws an ROC curve for the training set (for comparison)
#   - All results are saved to a folder named "validation_results"
# ====================================================================

# ------------------------- Load packages -------------------------
suppressPackageStartupMessages({
  library(pROC)
  library(rms)
})

# ------------------------- Path setup -------------------------
working_directory <- getwd()
setwd(working_directory)

# Output directory: validation results
out_dir <- "validation_results"
if (!dir.exists(out_dir)) dir.create(out_dir)

# Input files
train_expr_file <- "merged_after_batch_removal.csv"
gene_file       <- "Intersect_Genes.csv"
# External validation cohort (only GSE53757, as reported in the manuscript)
val_files <- c(
  GSE53757 = "GSE53757.csv"
)

# ------------------------- 读取基因列表 -------------------------
gene_tbl <- read.csv(gene_file, header = TRUE, check.names = FALSE,
                     stringsAsFactors = FALSE)
feature_genes <- trimws(as.character(gene_tbl[, 1]))
feature_genes <- feature_genes[feature_genes != ""]
cat(sprintf("读取到 %d 个诊断基因: %s\n",
            length(feature_genes), paste(feature_genes, collapse = ", ")))

# ------------------------- 读取训练集并重构模型 -------------------------
cat("\n==== 步骤1: 在训练集(merged_after_batch_removal)上重构模型 ====\n")

dat_train <- read.table(train_expr_file, sep = ",", header = TRUE,
                        row.names = 1, check.names = FALSE)
dat_train <- dat_train[gsub("-", "_", rownames(dat_train)), , drop = FALSE]

found_train <- feature_genes %in% rownames(dat_train)
if (any(!found_train)) {
  warning(paste0("以下基因在训练集中未找到: ",
                 paste(feature_genes[!found_train], collapse = ",")))
}
feature_genes <- feature_genes[found_train]

dat_train_filt <- dat_train[feature_genes, , drop = FALSE]
df_train <- as.data.frame(t(dat_train_filt))

sample_names_tr <- rownames(df_train)
groups_tr <- gsub("(.*)_([A-Za-z0-9]+)$", "\\2", sample_names_tr)
df_train$GroupType <- groups_tr

group_levels <- sort(unique(groups_tr))  # 字母顺序: con=0, tre=1
df_train$GroupType <- as.numeric(factor(df_train$GroupType,
                                         levels = group_levels)) - 1
cat(sprintf("训练集样本数: %d (Control=%d, Disease=%d)\n",
            nrow(df_train),
            sum(df_train$GroupType == 0), sum(df_train$GroupType == 1)))

ddist <- datadist(df_train)
options(datadist = "ddist")

model_vars <- setdiff(colnames(df_train), "GroupType")
reg_formula <- as.formula(paste("GroupType ~",
                                paste(model_vars, collapse = " + ")))

lrm_fit <- lrm(reg_formula, data = df_train, x = TRUE, y = TRUE)
cat("训练集模型构建完成。\n")

pred_train <- predict(lrm_fit, newdata = df_train, type = "fitted")

# ------------------------- 通用函数 -------------------------
# 计算AUC及其DeLong 95% CI
compute_auc_ci <- function(roc_obj) {
  auc_val <- as.numeric(auc(roc_obj))
  ci_res  <- ci.auc(roc_obj, conf.level = 0.95, method = "delong")
  ci_vec  <- as.numeric(ci_res)
  if (length(ci_vec) == 3) {
    CI_low  <- ci_vec[1]
    CI_high <- ci_vec[3]
  } else if (length(ci_vec) == 2) {
    CI_low  <- ci_vec[1]
    CI_high <- ci_vec[2]
  } else {
    stop("无法解析 ci.auc 返回值")
  }
  list(AUC = auc_val, CI_low = CI_low, CI_high = CI_high)
}

# 风格与原脚本一致：仅在图中显示 AUC（不显示 95% CI）
plot_roc_simple <- function(roc_obj, auc_info, main_title, outfile,
                            col_main = "red", lwd = 4) {
  pdf(outfile, width = 6, height = 6)
  par(mar = c(5, 6, 4, 2) + 0.1, cex = 1.3)

  plot(1 - roc_obj$specificities, roc_obj$sensitivities, type = "l",
       col = col_main, lwd = lwd, lty = 1,
       xlab = expression("1 - Specificity"),
       ylab = "Sensitivity",
       main = main_title,
       cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45,
       xlim = c(0, 1), ylim = c(0, 1))
  abline(0, 1, lty = 2, col = "gray70", lwd = 2)

  legend("bottomright",
         legend = sprintf("AUC = %.3f", auc_info$AUC),
         col = col_main, lwd = lwd, lty = 1, bty = "n", cex = 1.2)
  dev.off()
}

# ------------------------- 步骤2: 训练集 ROC -------------------------
cat("\n==== 步骤2: 训练集 ROC(原脚本未展示)====\n")

roc_train <- roc(df_train$GroupType, pred_train,
                 levels = c(0, 1), direction = "<")
auc_train <- compute_auc_ci(roc_train)

cat(sprintf("训练集 AUC = %.3f (95%% CI: %.3f - %.3f)\n",
            auc_train$AUC, auc_train$CI_low, auc_train$CI_high))

plot_roc_simple(
  roc_train, auc_train,
  main_title = "Model ROC Curve",
  outfile    = file.path(out_dir, "Training_ROC.pdf"),
  col_main   = "red"
)
cat(">> 训练集 ROC -> 验证结果/Training_ROC.pdf\n")

# ------------------------- 步骤3-N: 外部验证（统一循环）-------------------------
cat("\n==== 步骤3: 外部验证集 ROC ====\n")

# Colour scheme for each cohort
val_colors <- c(
  GSE53757 = "forestgreen"
)

# 用于汇总
cohort_names   <- c("Training Set", names(val_files))
n_samples      <- c(nrow(df_train), rep(NA, length(val_files)))
n_control      <- c(sum(df_train$GroupType == 0), rep(NA, length(val_files)))
n_disease      <- c(sum(df_train$GroupType == 1), rep(NA, length(val_files)))
auc_vals       <- c(auc_train$AUC, rep(NA, length(val_files)))
ci_low_vals    <- c(auc_train$CI_low, rep(NA, length(val_files)))
ci_high_vals   <- c(auc_train$CI_high, rep(NA, length(val_files)))

val_results <- list()  # 存储各验证集ROC对象、预测概率等

for (i in seq_along(val_files)) {
  cohort_name <- names(val_files)[i]
  val_file    <- val_files[i]
  cat(sprintf("\n---- %s ----\n", cohort_name))

  dat_v <- read.table(val_file, sep = ",", header = TRUE,
                      row.names = 1, check.names = FALSE)
  dat_v <- dat_v[gsub("-", "_", rownames(dat_v)), , drop = FALSE]

  found_v <- feature_genes %in% rownames(dat_v)
  if (any(!found_v)) {
    warning(paste0(cohort_name, " 中缺失基因: ",
                   paste(feature_genes[!found_v], collapse = ",")))
  }
  feature_v <- feature_genes[found_v]

  dat_v_filt <- dat_v[feature_v, , drop = FALSE]
  df_v <- as.data.frame(t(dat_v_filt))

  sample_names_v <- rownames(df_v)
  groups_v <- gsub("(.*)_([A-Za-z0-9]+)$", "\\2", sample_names_v)
  df_v$GroupType <- as.numeric(factor(groups_v, levels = group_levels)) - 1

  cat(sprintf("%s 样本数: %d (Control=%d, Disease=%d)\n",
              cohort_name, nrow(df_v),
              sum(df_v$GroupType == 0), sum(df_v$GroupType == 1)))

  pred_v <- predict(lrm_fit, newdata = df_v, type = "fitted")

  roc_v <- roc(df_v$GroupType, pred_v,
               levels = c(0, 1), direction = "<")
  auc_v <- compute_auc_ci(roc_v)
  cat(sprintf("%s AUC = %.3f (95%% CI: %.3f - %.3f)\n",
              cohort_name, auc_v$AUC, auc_v$CI_low, auc_v$CI_high))

  # 单独 ROC 曲线（图例只显示 AUC，无 95% CI）
  plot_roc_simple(
    roc_v, auc_v,
    main_title = sprintf("%s ROC Curve", cohort_name),
    outfile    = file.path(out_dir, sprintf("%s_ROC.pdf", cohort_name)),
    col_main   = val_colors[cohort_name]
  )
  cat(sprintf(">> %s ROC -> 验证结果/%s_ROC.pdf\n",
              cohort_name, cohort_name))

  # 保存预测概率
  pred_out <- data.frame(Sample = rownames(df_v),
                         GroupType = df_v$GroupType,
                         Predicted_Prob = pred_v)
  write.csv(pred_out,
            file = file.path(out_dir,
                             sprintf("%s_Predicted_Probabilities.csv",
                                     cohort_name)),
            row.names = FALSE)

  # 保存中间结果
  val_results[[cohort_name]] <- list(
    roc_obj = roc_v, auc = auc_v,
    n_samples = nrow(df_v),
    n_control = sum(df_v$GroupType == 0),
    n_disease = sum(df_v$GroupType == 1)
  )

  # 填入汇总
  idx <- i + 1
  n_samples[idx]    <- nrow(df_v)
  n_control[idx]    <- sum(df_v$GroupType == 0)
  n_disease[idx]    <- sum(df_v$GroupType == 1)
  auc_vals[idx]     <- auc_v$AUC
  ci_low_vals[idx]  <- auc_v$CI_low
  ci_high_vals[idx] <- auc_v$CI_high
}

# ------------------------- 步骤4: 多集合并 ROC 对比图 -------------------------
cat("\n==== 步骤4: 训练集 + 外部验证集 合并 ROC 对比图 ====\n")

pdf(file.path(out_dir, "Combined_ROC_All_Cohorts.pdf"),
    width = 8, height = 8)
par(mar = c(5, 6, 4, 2) + 0.1, cex = 1.3)

# 先画训练集作为基准
plot(1 - roc_train$specificities, roc_train$sensitivities, type = "l",
     col = "red", lwd = 3, lty = 1,
     xlab = expression("1 - Specificity"),
     ylab = "Sensitivity",
     main = "ROC Curves: Training & External Validations",
     cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.4,
     xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, lty = 2, col = "gray70", lwd = 2)

# 验证集叠加
legend_labels <- c(sprintf("Training (n=%d): AUC=%.3f",
                           nrow(df_train), auc_train$AUC))
legend_cols   <- c("red")

for (cohort_name in names(val_results)) {
  vr <- val_results[[cohort_name]]
  lines(1 - vr$roc_obj$specificities, vr$roc_obj$sensitivities,
        col = val_colors[cohort_name], lwd = 3, lty = 1)
  legend_labels <- c(legend_labels,
                     sprintf("%s (n=%d): AUC=%.3f",
                             cohort_name, vr$n_samples, vr$auc$AUC))
  legend_cols   <- c(legend_cols, val_colors[cohort_name])
}

legend("bottomright",
       legend = legend_labels,
       col    = legend_cols,
       lwd    = 3, lty = 1, bty = "n", cex = 1.05)
dev.off()
cat(">> 合并 ROC 对比图 -> 验证结果/Combined_ROC_All_Cohorts.pdf\n")

# ------------------------- 步骤5: 汇总 AUC 结果表 -------------------------
summary_df <- data.frame(
  Cohort      = cohort_names,
  N_Samples   = n_samples,
  N_Control   = n_control,
  N_Disease   = n_disease,
  AUC         = round(auc_vals, 4),
  CI_Lower_95 = round(ci_low_vals, 4),
  CI_Upper_95 = round(ci_high_vals, 4)
)
write.csv(summary_df,
          file = file.path(out_dir, "AUC_Summary_All_Cohorts.csv"),
          row.names = FALSE)

cat("\n==== AUC 汇总 ====\n")
print(summary_df, row.names = FALSE)
cat("\n==== 外部验证全流程执行完毕!====\n")
cat(">> 所有结果已保存至: 验证结果 文件夹\n")