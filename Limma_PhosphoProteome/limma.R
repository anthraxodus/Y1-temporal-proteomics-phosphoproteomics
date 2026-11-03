# =====================================================
# STEP 0: Load libraries
# =====================================================
library(limma)
library(openxlsx)

rm(list = ls())


# =====================================================
# STEP 1: Load and inspect data
# =====================================================
# Load the phospho dataset
phospho_clean <- read.delim(
  "Data/05_Phospho (STY)Sites_PIF_0.5_DownstreamReady.tsv",
  header = TRUE,
  check.names = FALSE,
  sep = "\t"
)

# Print first rows and column names to check
head(phospho_clean)
colnames(phospho_clean)  # all 72 intensity columns + 27 feature columns

# =====================================================
# STEP 2: Select intensity columns (samples)
# =====================================================
# First, get all columns that start with "Intensity "
intensity_cols <- grep("^Intensity ", colnames(phospho_clean))

# Exclude reference columns ("Intensity ref_")
ref_cols <- grep("^Intensity ref_", colnames(phospho_clean))
if(length(ref_cols) > 0){
  intensity_cols <- setdiff(intensity_cols, ref_cols)
  message("Removed ", length(ref_cols), " reference columns: ", 
          paste(colnames(phospho_clean)[ref_cols], collapse=", "))
}

# Convert to matrix
quant <- as.matrix(phospho_clean[, intensity_cols])
dim(quant)  # should be n_sites x (72 - number of ref columns)

# =====================================================
# STEP 3: Prepare feature annotations
# =====================================================
# Keep all original features for master table
feature_data <- phospho_clean[, 1:27]

# Assign unique IDs (from "T: Unique identifier")
rownames(quant) <- phospho_clean[,"T: Unique identifier"]

# =====================================================
# STEP 4: Prepare metadata for samples
# =====================================================
# Remove "Intensity " prefix
sample_names <- gsub("Intensity ", "", colnames(quant))
sample_split <- strsplit(sample_names, "_")

# Create metadata dataframe
metadata <- data.frame(
  condition = sapply(sample_split, "[[", 1),
  time = sapply(sample_split, "[[", 2),
  replicate = sapply(sample_split, "[[", 3),
  stringsAsFactors = FALSE
)
rownames(metadata) <- colnames(quant)
head(metadata)

# =====================================================
# STEP 5: Define groups for contrasts
# =====================================================
# Combine condition + time
grps <- paste0(metadata$condition, "_", metadata$time)

# Optional: filter grps, replicate selection, etc.
# grps_filtered <- grps

# =====================================================
# STEP 6: Design matrix for Limma
# =====================================================
design <- model.matrix(~0 + factor(grps))  # ~0 to remove intercept
colnames(design) <- levels(factor(grps))
head(design)

View(design)

# =====================================================
# STEP 7: Fit linear model
# =====================================================
fit <- lmFit(quant, design)
fit <- eBayes(fit)

# =====================================================
# STEP 8: Define contrasts manually
# =====================================================
contrast_list <- list(
  list(group1="FGF2minusFBS_5s",        group2="CtrlminusFBS_0m",  title="FGF2 - FBS 5 Seconds",       string="A."),
  list(group1="FGF2minusFBS_5m",        group2="CtrlminusFBS_5m",  title="FGF2 - FBS 5 Minutes",       string="B."),
  list(group1="FGF2minusFBS_10m",       group2="CtrlminusFBS_10m", title="FGF2 - FBS 10 Minutes",      string="C."),
  list(group1="FGF2minusFBS_30m",       group2="CtrlminusFBS_30m", title="FGF2 - FBS 30 Minutes",      string="D."),
  list(group1="FGF2plusFBS_5s",         group2="CtrlplusFBS_0m",   title="FGF2 + FBS 5 Seconds",       string="E."),
  list(group1="FGF2plusFBS_5m",         group2="CtrlplusFBS_5m",   title="FGF2 + FBS 5 Minutes",       string="F."),
  list(group1="FGF2plusFBS_10m",        group2="CtrlplusFBS_10m",  title="FGF2 + FBS 10 Minutes",      string="G."),
  list(group1="FGF2plusFBS_30m",        group2="CtrlplusFBS_30m",  title="FGF2 + FBS 30 Minutes",      string="H."),
  list(group1="SOSionFGF2minusFBS_5s",  group2="CtrlminusFBS_0m",  title="SOSi/FGF2 - FBS 5 Seconds",  string="I."),
  list(group1="SOSionFGF2minusFBS_5m",  group2="CtrlminusFBS_5m",  title="SOSi/FGF2 - FBS 5 Minutes",  string="J."),
  list(group1="SOSionFGF2minusFBS_10m", group2="CtrlminusFBS_10m", title="SOSi/FGF2 - FBS 10 Minutes", string="K."),
  list(group1="SOSionFGF2minusFBS_30m", group2="CtrlminusFBS_30m", title="SOSi/FGF2 - FBS 30 Minutes", string="L."),
  list(group1="SOSionFGF2plusFBS_5s",   group2="CtrlplusFBS_0m",   title="SOSi/FGF2 + FBS 5 Seconds",  string="M."),
  list(group1="SOSionFGF2plusFBS_5m",   group2="CtrlplusFBS_5m",   title="SOSi/FGF2 + FBS 5 Minutes",  string="N."),
  list(group1="SOSionFGF2plusFBS_10m",  group2="CtrlplusFBS_10m",  title="SOSi/FGF2 + FBS 10 Minutes", string="O."),
  list(group1="SOSionFGF2plusFBS_30m",  group2="CtrlplusFBS_30m",  title="SOSi/FGF2 + FBS 30 Minutes", string="P.")
)


View(grps)

# Prepare folder structure
general_dir <- "Outputs/limma_phospho/general"
dir.create(general_dir, recursive = TRUE, showWarnings = FALSE)

contrasts_base_dir <- "Outputs/limma_phospho/contrasts"
dir.create(contrasts_base_dir, recursive = TRUE, showWarnings = FALSE)

kinase_ready_cutoff <- file.path(contrasts_base_dir, "Z.Kinase-library-ready/cutoff")
kinase_ready_nocutoff <- file.path(contrasts_base_dir, "Z.Kinase-library-ready/no_cutoff")
dir.create(kinase_ready_cutoff, recursive = TRUE, showWarnings = FALSE)
dir.create(kinase_ready_nocutoff, recursive = TRUE, showWarnings = FALSE)

# Initialize master table with features + intensities
master_table <- cbind(feature_data, quant)

# =====================================================
# STEP 9: Loop through contrasts (corrected)
# =====================================================
# Definir thresholds
lfc_threshold <- log2(1.5)  # log2 fold change cutoff
pval_threshold <- 0.05      # adjusted p-value cutoff

# Loop pelos contrastes
for(contrast in contrast_list){
  
  # Extract contrast info
  g1 <- contrast$group1
  g2 <- contrast$group2
  title <- contrast$title
  string <- contrast$string
  
  message("Processing contrast: ", title)
  
  # Make safe folder name
  safe_title <- gsub("[/\\:]", "_", paste0(string, title))
  contrast_dir <- file.path(contrasts_base_dir, safe_title)
  dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create contrast vector
  contrast_str <- paste0(g1, "-", g2)
  cont_matrix <- makeContrasts(contrasts = contrast_str, levels=design)
  
  # Fit contrast and compute statistics
  fit2 <- contrasts.fit(fit, cont_matrix)
  fit2 <- eBayes(fit2)
  
  res <- topTable(fit2, coef=contrast_str, number=Inf)
  
  # Reorder rows to match master_table
  res <- res[match(rownames(master_table), rownames(res)), ]
  
  # Calculate -log10 adj P-value
  res$logP <- -log10(res$adj.P.Val)
  
  # Identify Up/Down and significant sites (Kinase library criteria)
  up <- which(res$adj.P.Val < pval_threshold & res$logFC > lfc_threshold)
  down <- which(res$adj.P.Val < pval_threshold & res$logFC < -lfc_threshold)
  UpDown <- rep("", nrow(res))
  UpDown[up] <- "Up"
  UpDown[down] <- "Down"
  significant <- ifelse(res$adj.P.Val < pval_threshold & abs(res$logFC) > lfc_threshold, "+", "")
  
  # Combine results with master table for this contrast
  contrast_table <- cbind(master_table[,1:27], res[, c("logFC","P.Value","adj.P.Val")])
  contrast_table$UpDown <- UpDown
  contrast_table$significant <- significant
  
  # Save contrast tables
  write.csv(contrast_table, file.path(contrast_dir, paste0(safe_title, ".csv")), row.names=FALSE)
  write.csv(contrast_table[contrast_table$significant=="+", ], 
            file.path(contrast_dir, paste0(safe_title, "_SignificantOnly.csv")), row.names=FALSE)
  openxlsx::write.xlsx(contrast_table[contrast_table$significant=="+", ],
                       file.path(contrast_dir, paste0(safe_title, "_SignificantOnly.xlsx")), rowNames=FALSE)
  
  # Save for Kinase-library-ready (with cutoff)
  write.csv(contrast_table[contrast_table$significant=="+", ], 
            file.path(kinase_ready_cutoff, paste0(safe_title, ".csv")), row.names=FALSE)
  
  # Save for Kinase-library-ready (without cutoff = all)
  write.csv(contrast_table, 
            file.path(kinase_ready_nocutoff, paste0(safe_title, ".csv")), row.names=FALSE)
  
  # Volcano plot
  png(filename = file.path(contrast_dir, paste0(safe_title, ".png")))
  plot(res$logFC, res$logP, main=title, xlab="Log2 Fold Change", ylab="-log10(adj P-value)",
       pch=20, col="grey")
  points(res$logFC[up], res$logP[up], col="red", pch=16)
  points(res$logFC[down], res$logP[down], col="blue", pch=16)
  text(x=min(res$logFC), y=max(res$logP)*0.95, labels=paste0("Up: ", length(up)), col="red", adj=0)
  text(x=min(res$logFC), y=max(res$logP)*0.9, labels=paste0("Down: ", length(down)), col="blue", adj=0)
  dev.off()
  
  # Add results to master_table (general)
  master_table[[paste0("logFC_", string)]] <- res$logFC
  master_table[[paste0("P.Value_", string)]] <- res$P.Value
  master_table[[paste0("adj.P.Val_", string)]] <- res$adj.P.Val
  master_table[[paste0("significant_", string)]] <- significant
}

# =====================================================
# STEP 10: General significance column
# =====================================================
sig_cols <- grep("^significant_", colnames(master_table), value=TRUE)
master_table$significant_any <- apply(master_table[, sig_cols], 1, function(x){
  if(any(x == "+")) return("+") else return("")
})

# =====================================================
# STEP 11: Save master table (general)
# =====================================================
write.csv(master_table, file.path(general_dir, "06_Phospho (STY)Sites_PIF_0.5_limmaSignificant.csv"), row.names=FALSE)
openxlsx::write.xlsx(master_table, file.path(general_dir, "06_Phospho (STY)Sites_PIF_0.5_limmaSignificant.xlsx"), rowNames=FALSE)

                                           