options(timeout = 1800)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  message("...")
} else {
  BiocManager::install("TCGAbiolinks", update = FALSE, ask = FALSE)
}

library(TCGAbiolinks)
library(SummarizedExperiment)

tcga_projects <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD",
  "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC",
  "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
  "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ",
  "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)

dir.create("tcga_counts", showWarnings = FALSE)

for (project in tcga_projects) {
  message("Processing: ", project)

  tryCatch({
    query <- GDCquery(
      project = project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    GDCdownload(query, method = "api", files.per.chunk = 20)
    data <- GDCprepare(query)

    sample_type <- data$sample_type
    group_label <- ifelse(sample_type == "Solid Tissue Normal", "NORMAL", "TUMOR")

    if (sum(group_label == "NORMAL") < 3) {
      message(project, ": fewer than 3 normal samples")
    }

    counts <- assay(data, "unstranded")
    colnames(counts) <- paste0(colnames(counts), "_", group_label)

    counts_df <- as.data.frame(counts)
    counts_df <- cbind(gene_id = rownames(counts_df), counts_df)

    out_file <- file.path(
      "tcga_counts",
      paste0(gsub("-", "_", project), "_counts.csv")
    )
    write.csv(counts_df, file = out_file, row.names = FALSE)
    message(project, ": wrote ", out_file, " (", ncol(counts), " samples)")

  },
  error = function(e) {
    message(project, ": failed - ", conditionMessage(e))
  })
}