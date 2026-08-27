# ========================================================================
# Required Libraries
# ========================================================================
library(GOplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(cowplot)

# ==============================================================================
# 18_go_kegg_enrichment.R
# Functional enrichment analysis (GO and KEGG) for hub candidate genes.
# ==============================================================================

# ========================================================================
# Set working directory (relative path; no absolute paths used)
# ========================================================================
workDir <- getwd()
setwd(workDir)

# ========================================================================
# Auto-detect and load gene files
# ========================================================================
message("Scanning for gene files starting with '05_' and ending with '.txt'...")

# List all files matching the pattern
gene_files <- list.files(pattern = "^05_.*\\.txt$", full.names = FALSE)

if (length(gene_files) == 0) {
  stop("No gene files found matching pattern '05_*.txt'")
}

message(sprintf("Found %d gene file(s):", length(gene_files)))
for (i in seq_along(gene_files)) {
  message(sprintf("  [%d] %s", i, gene_files[i]))
}

# Use the first matching file
input_gene_file <- gene_files[1]
message(sprintf("\nUsing gene file: %s", input_gene_file))

# ========================================================================
# Load data and ID conversion
# ========================================================================
message("Loading gene data...")

# Read the gene list from txt file
gene_symbols <- read.table(input_gene_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
gene_symbols <- trimws(gene_symbols)
gene_symbols <- gene_symbols[gene_symbols != ""]
gene_symbols <- unique(gene_symbols)

message(sprintf("Total genes: %d", length(gene_symbols)))

# ========================================================================
# Check for logFC file (optional)
# ========================================================================
message("Checking for differential expression file...")

# Create a mapping of gene to logFC values
gene_logFC_map <- data.frame(Gene = character(), logFC = numeric(), stringsAsFactors = FALSE)
has_logFC <- FALSE

# Try to find DE file
de_files <- c("DE_significant_genes.csv", "differential_expression.csv", "DE_results.csv")
de_file_found <- NULL

for (de_file in de_files) {
  if (file.exists(de_file)) {
    de_file_found <- de_file
    break
  }
}

if (!is.null(de_file_found)) {
  message(sprintf("Found differential expression file: %s", de_file_found))
  tryCatch({
    de_data <- read.csv(de_file_found, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

    message(sprintf("File loaded successfully with %d rows and %d columns", nrow(de_data), ncol(de_data)))
    message(sprintf("Column names: %s", paste(colnames(de_data), collapse = ", ")))

    # Check if Gene and logFC columns exist
    if ("Gene" %in% colnames(de_data) && "logFC" %in% colnames(de_data)) {
      # Extract Gene and logFC columns
      gene_logFC_map <- de_data[, c("Gene", "logFC")]
      gene_logFC_map <- gene_logFC_map[!is.na(gene_logFC_map$Gene) & !is.na(gene_logFC_map$logFC), ]
      message(sprintf("Extracted %d genes with logFC values", nrow(gene_logFC_map)))
      has_logFC <- TRUE
    } else {
      message(sprintf("Warning: Missing Gene or logFC column in %s", de_file_found))
      message(sprintf("Available columns: %s", paste(colnames(de_data), collapse = ", ")))
    }
  }, error = function(e) {
    message(sprintf("Error reading %s: %s", de_file_found, e$message))
  })
} else {
  message("No differential expression file found. Will use default logFC = 1 for all genes.")
}

# Remove duplicates, keeping the first occurrence
if (has_logFC) {
  gene_logFC_map <- gene_logFC_map[!duplicated(gene_logFC_map$Gene), ]
  message(sprintf("Total unique genes with logFC values: %d", nrow(gene_logFC_map)))
}

# Convert gene symbols to Entrez IDs
message("Converting gene symbols to Entrez IDs...")
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

message(sprintf("Successfully converted: %d genes", nrow(entrez_ids)))

# Merge with logFC values
message("Merging logFC values with gene list...")
logFC_values <- rep(1, nrow(entrez_ids))
names(logFC_values) <- entrez_ids$ENTREZID

# Update logFC values from the diff files
for (i in 1:nrow(entrez_ids)) {
  gene_symbol <- entrez_ids$SYMBOL[i]
  entrez_id <- entrez_ids$ENTREZID[i]

  # Look up logFC value
  logFC_idx <- which(gene_logFC_map$Gene == gene_symbol)
  if (length(logFC_idx) > 0) {
    logFC_values[as.character(entrez_id)] <- gene_logFC_map$logFC[logFC_idx[1]]
  }
}

# Create sorted gene list (descending by logFC)
genelist <- sort(logFC_values, decreasing = TRUE)

message(sprintf("Gene list prepared with %d genes", length(genelist)))

# ========================================================================
# Create gene list for GOplot with actual logFC values
# ========================================================================
genelist_df <- data.frame(
  ID = entrez_ids$SYMBOL,
  logFC = rep(1, nrow(entrez_ids)),
  stringsAsFactors = FALSE
)

# Update logFC values in genelist_df
for (i in 1:nrow(genelist_df)) {
  gene_symbol <- genelist_df$ID[i]
  logFC_idx <- which(gene_logFC_map$Gene == gene_symbol)
  if (length(logFC_idx) > 0) {
    genelist_df$logFC[i] <- gene_logFC_map$logFC[logFC_idx[1]]
  }
}

row.names(genelist_df) <- genelist_df$ID

# ========================================================================
# Create gene color mapping based on logFC
# ========================================================================
message("Creating gene color mapping...")

# Function to get color based on logFC
get_gene_color <- function(logFC_value) {
  if (logFC_value > 0) {
    return("#FF6347")  # Orange-red for upregulated
  } else if (logFC_value < 0) {
    return("#6A5ACD")  # Blue-purple for downregulated
  } else {
    return("#808080")  # Gray for neutral
  }
}

# Create color vector for all genes
gene_colors <- sapply(genelist_df$logFC, get_gene_color)
names(gene_colors) <- genelist_df$ID

message(sprintf("Upregulated genes (orange-red): %d", sum(genelist_df$logFC > 0)))
message(sprintf("Downregulated genes (blue-purple): %d", sum(genelist_df$logFC < 0)))

# ========================================================================
# GO Enrichment Analysis
# ========================================================================
message("Performing GO enrichment analysis...")

ego_bp <- enrichGO(
  gene = names(genelist),
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ego_cc <- enrichGO(
  gene = names(genelist),
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ego_mf <- enrichGO(
  gene = names(genelist),
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pvalueCutoff = 0.05,
  readable = TRUE
)

message(sprintf("GO BP terms found: %d", nrow(ego_bp@result)))
message(sprintf("GO CC terms found: %d", nrow(ego_cc@result)))
message(sprintf("GO MF terms found: %d", nrow(ego_mf@result)))

# Combine all GO results - only include non-empty results
ego_list <- list()
if (nrow(ego_bp@result) > 0) {
  ego_bp_df <- ego_bp@result
  ego_bp_df$ONTOLOGY <- "BP"
  ego_list[[length(ego_list) + 1]] <- ego_bp_df
}

if (nrow(ego_cc@result) > 0) {
  ego_cc_df <- ego_cc@result
  ego_cc_df$ONTOLOGY <- "CC"
  ego_list[[length(ego_list) + 1]] <- ego_cc_df
}

if (nrow(ego_mf@result) > 0) {
  ego_mf_df <- ego_mf@result
  ego_mf_df$ONTOLOGY <- "MF"
  ego_list[[length(ego_list) + 1]] <- ego_mf_df
}

if (length(ego_list) > 0) {
  ego <- do.call(rbind, ego_list)
  row.names(ego) <- NULL
} else {
  message("Warning: No GO terms found!")
  ego <- data.frame(
    ONTOLOGY = character(0),
    ID = character(0),
    Description = character(0),
    geneID = character(0),
    p.adjust = numeric(0)
  )
}

# ========================================================================
# KEGG Enrichment Analysis
# ========================================================================
message("Performing KEGG pathway analysis...")

kegg <- enrichKEGG(
  gene = names(genelist),
  organism = "hsa",
  pvalueCutoff = 1
)

# ========================================================================
# Prepare data for visualization
# ========================================================================
message("Preparing visualization data...")

# Function to create circular data for GOplot with gene colors
create_goplot_circular_data <- function(go_result, genelist, go_type, gene_colors) {
  tryCatch({
    go_df <- data.frame(
      Category = go_type,
      ID = go_result$ID,
      Term = go_result$Description,
      Genes = gsub("/", ", ", go_result$geneID),
      adj_pval = go_result$p.adjust,
      stringsAsFactors = FALSE
    )

    if (nrow(go_df) > 0) {
      # Try to extract genes that exist in genelist
      genes_in_go <- unique(unlist(strsplit(gsub(", ", ",", go_df$Genes), ",")))
      genelist_subset <- genelist[genelist$ID %in% genes_in_go, ]

      if (nrow(genelist_subset) > 0) {
        # Add color information to genelist_subset
        genelist_subset$color <- gene_colors[genelist_subset$ID]

        circ <- circle_dat(go_df, genelist_subset)
        return(list(go = go_df, circ = circ, success = TRUE, colors = gene_colors[genes_in_go]))
      }
    }
    return(list(success = FALSE))
  }, error = function(e) {
    message(sprintf("Error preparing %s data: %s", go_type, e$message))
    return(list(success = FALSE))
  })
}

# Prepare GO BP data
message("Preparing GO BP (Biological Process) data...")
go_bp_data <- create_goplot_circular_data(ego_bp@result, genelist_df, "BP", gene_colors)

# Prepare GO CC data
message("Preparing GO CC (Cellular Component) data...")
go_cc_data <- create_goplot_circular_data(ego_cc@result, genelist_df, "CC", gene_colors)

# Prepare GO MF data
message("Preparing GO MF (Molecular Function) data...")
go_mf_data <- create_goplot_circular_data(ego_mf@result, genelist_df, "MF", gene_colors)

# Prepare KEGG data
message("Preparing KEGG pathway data...")
kegg_data <- list(success = FALSE)
if (!is.null(kegg) && nrow(kegg@result) > 0) {
  tryCatch({
    # Use p.adjust (adjusted p-value) for filtering
    kegg_filtered <- kegg@result[kegg@result$p.adjust < 0.5, ]

    if (nrow(kegg_filtered) > 0) {
      # Create gene ID mapping BEFORE creating kegg_df
      gene_symbol_mapping <- setNames(entrez_ids$SYMBOL, entrez_ids$ENTREZID)

      # Convert Entrez IDs in geneID column to gene symbols
      kegg_filtered_copy <- kegg_filtered
      kegg_filtered_copy$geneID <- sapply(kegg_filtered$geneID, function(entrez_str) {
        entrez_ids_vec <- unlist(strsplit(entrez_str, "/"))
        symbols <- gene_symbol_mapping[entrez_ids_vec]
        symbols <- symbols[!is.na(symbols)]
        paste(symbols, collapse = "/")
      })

      kegg_df <- data.frame(
        Category = "KEGG",
        ID = kegg_filtered_copy$ID,
        Term = kegg_filtered_copy$Description,
        Genes = gsub("/", ", ", kegg_filtered_copy$geneID),
        adj_pval = kegg_filtered_copy$p.adjust,
        stringsAsFactors = FALSE
      )

      if (nrow(kegg_df) > 0) {
        message(sprintf("Preparing KEGG data with %d pathways", nrow(kegg_df)))

        # Now simply use circle_dat with genelist_df (which has gene symbols)
        tryCatch({
          # Add color information to genelist_df for KEGG
          genelist_df_with_color <- genelist_df
          genelist_df_with_color$color <- gene_colors[genelist_df_with_color$ID]

          kegg_circ <- circle_dat(kegg_df, genelist_df_with_color)
          kegg_data <- list(go = kegg_df, circ = kegg_circ, success = TRUE, colors = gene_colors)
          message("KEGG circular data prepared successfully")
        }, error = function(e) {
          message(sprintf("Error in circle_dat for KEGG: %s", e$message))

          # Alternative: create subset matching approach
          genes_in_kegg <- unique(unlist(strsplit(gsub(", ", ",", kegg_df$Genes), ",")))
          genelist_subset <- genelist_df[genelist_df$ID %in% genes_in_kegg, ]

          if (nrow(genelist_subset) > 0) {
            message(sprintf("Using subset approach: %d genes matched", nrow(genelist_subset)))
            genelist_subset$color <- gene_colors[genelist_subset$ID]
            kegg_circ <- circle_dat(kegg_df, genelist_subset)
            kegg_data <- list(go = kegg_df, circ = kegg_circ, success = TRUE, colors = gene_colors[genes_in_kegg])
            message("KEGG circular data prepared successfully (subset method)")
          }
        })
      }
    } else {
      message(sprintf("No KEGG pathways with p.adjust < 0.5. Total pathways: %d", nrow(kegg@result)))
    }
  }, error = function(e) {
    message(sprintf("Error preparing KEGG data: %s", e$message))
  })
}

# ========================================================================
# Visualization
# ========================================================================
message("Generating visualizations...")

# Function to create GOCircle plot with gene colors
create_gocircle_plot <- function(data, filename, plot_type) {
  if (data$success && !is.null(data$circ)) {
    message(sprintf("Creating %s GOCircle plot...", plot_type))
    tryCatch({
      pdf(file = filename, width = 11, height = 6)
      GOCircle(data$circ, rad1 = 2.5, rad2 = 3.5, label.size = 4, nsub = 10)
      dev.off()
      message(sprintf("Saved: %s", filename))
      return(TRUE)
    }, error = function(e) {
      if (dev.cur() != 1) dev.off()
      message(sprintf("Error creating %s GOCircle plot: %s", plot_type, e$message))
      return(FALSE)
    })
  } else {
    message(sprintf("Skipping %s GOCircle plot (no results)", plot_type))
    return(FALSE)
  }
}

# Generate GOCircle plots for each GO type
create_gocircle_plot(go_bp_data, "GOCircle_BP.pdf", "BP")
create_gocircle_plot(go_cc_data, "GOCircle_CC.pdf", "CC")
create_gocircle_plot(go_mf_data, "GOCircle_MF.pdf", "MF")

# Generate GOCircle plot for KEGG
create_gocircle_plot(kegg_data, "GOCircle_KEGG.pdf", "KEGG")

# ========================================================================
# Generate Bar plots and Bubble plots
# ========================================================================
message("Generating bar plots and bubble plots...")

# Create function to generate combined bar and bubble plots
create_go_plots <- function(go_list, plot_type_name) {
  if (length(go_list) == 0) {
    message(sprintf("No GO results to plot for %s", plot_type_name))
    return()
  }

  # Force specific output file names
  barplot_file <- sprintf("GO_%s_barplot.pdf", plot_type_name)
  bubble_file <- sprintf("GO_%s_bubble.pdf", plot_type_name)

  # Create bar plot
  tryCatch({
    message(sprintf("Creating %s bar plot...", plot_type_name))

    # Create plots for each ontology and combine
    plot_list <- list()

    # Generate plots for BP, CC, MF in order, showing top 10 by p-value
    for (ont in c("BP", "CC", "MF")) {
      if (ont %in% names(go_list) && nrow(go_list[[ont]]@result) > 0) {
        # Get original data and sort
        original_df <- go_list[[ont]]@result
        sorted_df <- original_df[order(original_df$pvalue), ]
        n_terms <- min(10, nrow(sorted_df))
        top10_df <- head(sorted_df, n_terms)

        message(sprintf("  - Creating barplot for %s (showing %d out of %d total terms)",
                       ont, n_terms, nrow(original_df)))
        message(sprintf("    First term: %s (pvalue: %.2e)", top10_df$Description[1], top10_df$pvalue[1]))
        message(sprintf("    Last term: %s (pvalue: %.2e)", top10_df$Description[n_terms], top10_df$pvalue[n_terms]))

        # Create ggplot bar plot
        top10_df$Description <- factor(top10_df$Description, levels = rev(top10_df$Description))

        # Calculate negative log p-value for color gradient
        top10_df$neg_log_pvalue <- -log10(top10_df$pvalue)

        p <- ggplot(top10_df, aes(x = Count, y = Description, fill = neg_log_pvalue)) +
          geom_bar(stat = "identity") +
          scale_fill_gradientn(colors = c("#8B4789", "#2E8B9E", "#2FAA4F", "#FDB462", "#E74C3C")) +
          theme_bw() +
          theme(axis.text.y = element_text(size = 10)) +
          labs(x = "Gene Count", y = "", title = paste0(ont, " (Top ", n_terms, " by p-value)"),
               fill = "-log10(p-value)")

        plot_list[[ont]] <- p
      }
    }

    message(sprintf("Total plots created for barplot: %d", length(plot_list)))

    # Adjust PDF height based on number of plots
    pdf_height <- max(8, length(plot_list) * 6)
    pdf(file = barplot_file, width = 12, height = pdf_height)

    if (length(plot_list) > 0) {
      go_bar <- cowplot::plot_grid(plotlist = plot_list, ncol = 1, align = "v", axis = "lr", rel_heights = rep(1, length(plot_list)))
      print(go_bar)
    }
    dev.off()
    message(sprintf("Saved: %s (contains %d categories)", barplot_file, length(plot_list)))
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    message(sprintf("Error creating %s bar plot: %s", plot_type_name, e$message))
    traceback()
  })

  # Create bubble plot
  tryCatch({
    message(sprintf("Creating %s bubble plot...", plot_type_name))

    plot_list <- list()

    # Generate plots for BP, CC, MF in order, showing top 10 by p-value
    for (ont in c("BP", "CC", "MF")) {
      if (ont %in% names(go_list) && nrow(go_list[[ont]]@result) > 0) {
        # Get original data and sort
        original_df <- go_list[[ont]]@result
        sorted_df <- original_df[order(original_df$pvalue), ]
        n_terms <- min(10, nrow(sorted_df))
        top10_df <- head(sorted_df, n_terms)

        message(sprintf("  - Creating bubble plot for %s (showing %d out of %d total terms)",
                       ont, n_terms, nrow(original_df)))
        message(sprintf("    First term: %s (pvalue: %.2e)", top10_df$Description[1], top10_df$pvalue[1]))
        message(sprintf("    Last term: %s (pvalue: %.2e)", top10_df$Description[n_terms], top10_df$pvalue[n_terms]))

        # Create ggplot bubble plot
        top10_df$Description <- factor(top10_df$Description, levels = rev(top10_df$Description))

        # Calculate GeneRatio numeric value
        top10_df$GeneRatio_num <- sapply(strsplit(as.character(top10_df$GeneRatio), "/"), function(x) {
          as.numeric(x[1]) / as.numeric(x[2])
        })

        # Calculate negative log p-value for color gradient
        top10_df$neg_log_pvalue <- -log10(top10_df$pvalue)

        p <- ggplot(top10_df, aes(x = GeneRatio_num, y = Description)) +
          geom_point(aes(size = Count, color = neg_log_pvalue)) +
          scale_color_gradientn(colors = c("#8B4789", "#2E8B9E", "#2FAA4F", "#FDB462", "#E74C3C")) +
          scale_size_continuous(range = c(3, 8)) +
          theme_bw() +
          theme(axis.text.y = element_text(size = 10)) +
          labs(x = "Gene Ratio", y = "", title = paste0(ont, " (Top ", n_terms, " by p-value)"),
               color = "-log10(p-value)", size = "Count")

        plot_list[[ont]] <- p
      }
    }

    message(sprintf("Total plots created for bubble: %d", length(plot_list)))

    # Adjust PDF height based on number of plots
    pdf_height <- max(8, length(plot_list) * 6)
    pdf(file = bubble_file, width = 12, height = pdf_height)

    if (length(plot_list) > 0) {
      go_bubble <- cowplot::plot_grid(plotlist = plot_list, ncol = 1, align = "v", axis = "lr", rel_heights = rep(1, length(plot_list)))
      print(go_bubble)
    }
    dev.off()
    message(sprintf("Saved: %s (contains %d categories)", bubble_file, length(plot_list)))
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    message(sprintf("Error creating %s bubble plot: %s", plot_type_name, e$message))
    traceback()
  })
}

# Prepare GO list for plotting
go_list <- list()
if (nrow(ego_bp@result) > 0) {
  go_list[["BP"]] <- ego_bp
  message(sprintf("Adding BP to go_list: %d terms", nrow(ego_bp@result)))
}
if (nrow(ego_cc@result) > 0) {
  go_list[["CC"]] <- ego_cc
  message(sprintf("Adding CC to go_list: %d terms", nrow(ego_cc@result)))
}
if (nrow(ego_mf@result) > 0) {
  go_list[["MF"]] <- ego_mf
  message(sprintf("Adding MF to go_list: %d terms", nrow(ego_mf@result)))
}

message(sprintf("Total categories in go_list: %d", length(go_list)))
message(sprintf("Categories: %s", paste(names(go_list), collapse=", ")))

# Generate combined GO plots
if (length(go_list) > 0) {
  create_go_plots(go_list, "Combined")
} else {
  message("WARNING: No GO results available for plotting!")
}

# KEGG Bar plot
if (!is.null(kegg) && nrow(kegg@result) > 0) {
  tryCatch({
    message("Creating KEGG bar plot...")
    pdf(file = "KEGG_barplot.pdf", width = 10, height = 8)
    kegg_bar <- barplot(
      kegg,
      drop = TRUE,
      showCategory = 15,
      label_format = 50,
      color = "p.adjust"
    ) +
      scale_fill_gradientn(colors = c("#8B4789", "#2E8B9E", "#2FAA4F", "#FDB462", "#E74C3C"))
    print(kegg_bar)
    dev.off()
    message("Saved: KEGG_barplot.pdf")
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    message(sprintf("Error creating KEGG bar plot: %s", e$message))
  })
}

# KEGG Bubble plot
if (!is.null(kegg) && nrow(kegg@result) > 0) {
  tryCatch({
    message("Creating KEGG bubble plot...")
    pdf(file = "KEGG_bubble.pdf", width = 10, height = 8)
    kegg_bubble <- dotplot(
      kegg,
      showCategory = 15,
      orderBy = "GeneRatio",
      label_format = 50,
      color = "p.adjust"
    ) +
      scale_color_gradientn(colors = c("#E74C3C", "#FDB462", "#2FAA4F", "#2E8B9E", "#8B4789"))
    print(kegg_bubble)
    dev.off()
    message("Saved: KEGG_bubble.pdf")
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    message(sprintf("Error creating KEGG bubble plot: %s", e$message))
  })
}

# ========================================================================
# Save enrichment results
# ========================================================================
message("Saving enrichment analysis results...")

# Save GO results
if (nrow(ego_bp@result) > 0) {
  write.csv(ego_bp@result, file = "GO_BP_Results.csv", row.names = FALSE)
  message("Saved: GO_BP_Results.csv")
}

if (nrow(ego_cc@result) > 0) {
  write.csv(ego_cc@result, file = "GO_CC_Results.csv", row.names = FALSE)
  message("Saved: GO_CC_Results.csv")
}

if (nrow(ego_mf@result) > 0) {
  write.csv(ego_mf@result, file = "GO_MF_Results.csv", row.names = FALSE)
  message("Saved: GO_MF_Results.csv")
}

# Save KEGG results
if (!is.null(kegg) && nrow(kegg@result) > 0) {
  write.csv(kegg@result, file = "KEGG_Results.csv", row.names = FALSE)
  message("Saved: KEGG_Results.csv")
}

# Save converted gene ID list
write.csv(entrez_ids, file = "Converted_Gene_IDs.csv", row.names = FALSE)
message("Saved: Converted_Gene_IDs.csv")

# Save gene logFC mapping
write.csv(genelist_df, file = "Gene_logFC_Mapping.csv", row.names = FALSE)
message("Saved: Gene_logFC_Mapping.csv")

# ========================================================================
# Summary
# ========================================================================
message("\n========== ANALYSIS SUMMARY ==========")
message(sprintf("Input gene file: %s", input_gene_file))
message(sprintf("Total input genes: %d", length(gene_symbols)))
message(sprintf("Successfully converted to Entrez IDs: %d", nrow(entrez_ids)))
if (has_logFC) {
  message(sprintf("Genes with logFC values: %d", nrow(gene_logFC_map)))
  message("\nGene Regulation Summary:")
  message(sprintf("  - Upregulated genes (logFC > 0): %d", sum(genelist_df$logFC > 0)))
  message(sprintf("  - Downregulated genes (logFC < 0): %d", sum(genelist_df$logFC < 0)))
} else {
  message("No logFC data available - all genes treated equally")
}
message("\nGO Enrichment Results:")
message(sprintf("  - BP (Biological Process): %d terms", nrow(ego_bp@result)))
message(sprintf("  - CC (Cellular Component): %d terms", nrow(ego_cc@result)))
message(sprintf("  - MF (Molecular Function): %d terms", nrow(ego_mf@result)))

if (!is.null(kegg)) {
  message("\nKEGG Pathway Results:")
  message(sprintf("  - Pathways: %d", nrow(kegg@result)))
}

message("\nAll analysis complete!")
message(sprintf("Results saved in: %s", getwd()))
