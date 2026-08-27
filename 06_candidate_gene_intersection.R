# ==============================================================================
# 06_candidate_gene_intersection.R
# Take the intersection of candidate gene lists produced by SVM-RFE
# (04_svm_rfe_feature_selection.R) and Random Forest
# (05_random_forest_gene_selection.R). The overlapping genes are the
# candidate core genes reported in the manuscript (Methods §2).
# The script also generates a Venn diagram.
# ==============================================================================

# --- Load necessary packages ---
library(ggvenn)
library(gridExtra)
library(ggplot2)

# --- 定义字符串重复函数 ---
`%R%` <- function(x, n) {
  paste(rep(x, n), collapse = "")
}

# --- Set working directory (relative; no absolute paths used) ---
setwd(getwd())

# --- 创建输出文件夹 ---
outputDir <- "output_folder"
if (!dir.exists(outputDir)) dir.create(outputDir)

# --- 获取所有txt和csv文件（不包含disease.txt和disease.csv）---
txt_files <- list.files(pattern="\\.txt$")
csv_files <- list.files(pattern="\\.csv$")

# 合并两种格式的文件
files <- c(txt_files, csv_files)
files <- files[files != "disease.txt" & files != "disease.csv"]

cat("=" %R% 60, "\n")
cat("检测到的输入文件\n")
cat("=" %R% 60, "\n")
cat(sprintf("TXT文件: %d 个\n", length(txt_files)))
cat(sprintf("CSV文件: %d 个\n", length(csv_files)))
cat(sprintf("总文件数: %d 个\n\n", length(files)))

# --- 读取所有文件中的基因信息，整理到geneList ---
geneList <- list()
for (inputFile in files) {
  # 根据文件扩展名选择读取方式
  if (tools::file_ext(inputFile) == "txt") {
    # TXT文件：直接读取每一行
    geneNames <- readLines(inputFile, warn=FALSE)
    geneNames <- trimws(geneNames)         # 去掉前后空格
    geneNames <- geneNames[geneNames != ""]# 去掉空字符串和空行

  } else if (tools::file_ext(inputFile) == "csv") {
    # CSV文件：读取为数据框，提取所有非空值
    data <- read.csv(inputFile, header = TRUE, check.names = FALSE)
    # 将所有列合并为一个向量，移除NA和空值
    geneNames <- as.character(unlist(data))
    geneNames <- trimws(geneNames)
    geneNames <- geneNames[geneNames != "" & !is.na(geneNames)]
  }

  # 去除重复
  uniqGene <- unique(geneNames)
  setName <- tools::file_path_sans_ext(basename(inputFile))  # 文件名去扩展名作为集合名
  geneList[[setName]] <- uniqGene

  cat(sprintf("✓ 已读取 %s，基因数: %d 个\n", inputFile, length(uniqGene)))
}

cat("\n")

# --- 计算所有基因名的并集 ---
unionGenes <- Reduce(union, geneList)
unionCount <- length(unionGenes)

cat(sprintf("并集基因总数: %d 个\n\n", unionCount))

# --- 绘制Venn图并保存 ---
pdf(file = file.path(outputDir, "venn.pdf"), width = 6, height = 6)
venn_plot <- ggvenn(
  geneList,
  show_percentage = TRUE,
  stroke_color = "white",
  stroke_size = 0.5,
  fill_color = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(geneList)], # 适应多个集合
  set_name_color = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(geneList)],
  set_name_size = 6,
  text_size = 4.5
)
# 添加并集个数
print(venn_plot)
dev.off()

setSizes <- sapply(geneList, length)
pdf(file = file.path(outputDir, "set_barplot.pdf"), width=6, height=4)
ggplot(data.frame(Set=names(setSizes), Size=setSizes), aes(x=Set, y=Size, fill=Set)) +
  geom_bar(stat="identity") +
  theme_minimal() + ylab("Gene Count") + xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

# --- 计算交集基因并输出为CSV ---
intersectGenes <- Reduce(intersect, geneList)
cat(sprintf("交集基因总数: %d 个\n", length(intersectGenes)))

write.csv(data.frame(Gene = intersectGenes),
          file = file.path(outputDir, "Intersect_Genes.csv"),
          row.names = FALSE)
cat(sprintf("交集基因已保存至: %s\n", file.path(outputDir, "Intersect_Genes.csv")))

cat("=" %R% 60, "\n")
cat("处理完成！已输出Venn图。\n")
cat("=" %R% 60, "\n")
