# Load required libraries
library(MSstats)
library(ggplot2)
library(reshape2)
library(stringr)
library(gplots)

# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital//Data/JohnsonLab/AndrewKurland/PolyIC-Timecourse/20230227-AK10-analysis/20230227-AK10-ph/msstats/")

# Set prefix for naming output files
prefix <- "20230321-AK10-ph"

# Read in modacc formatted Spectronaut output
raw <- read.delim(file = "../../20210907-AK10-ph/20210903-AK10p_Report-modacc.txt", quote = "")

# Convert to MSstats format
quant <- SpectronauttoMSstatsFormat(input = raw, intensity = "NormalizedPeakArea")


# Remove 8h samples and 24h_1

# Make a column combining Condition and Bio Replicate
quant$Condition_BioReplicate <- paste(quant$Condition, quant$BioReplicate)

# Make a list of samples to exclude
exclude <- c("x08hr 1", "x08hr 2", "x08hr 3", "x24hr 1")

# Filter out samples to exclude
quant <- quant[which(!(quant$Condition_BioReplicate %in% exclude)),]

# Make a list of conditions to include
include <- c("Control", "x05min", "x30min", "x02hr", "x24hr")

# Filter for include list
quant <- quant[which(quant$Condition %in% include),]

# Add a column with log2intensity
quant$Log2Intensity <- log(quant$Intensity) / log(2)

# Calculate median log2intensity by run
median_log2intensity <- aggregate(Log2Intensity ~ Run, quant, median)

# Get the max median log2 intensity
max_median_log2intensity <- max(median_log2intensity$Log2Intensity)

# Calculate the distance from the max for each run
median_log2intensity$Diff <- max_median_log2intensity - median_log2intensity$Log2Intensity

# Add diff value to each run
for (i in 1:nrow(median_log2intensity)) {

  # Get current run and diff
  current_run <- median_log2intensity[i, "Run"]
  current_diff <- median_log2intensity[i, "Diff"]

  # Add diff
  quant[which(quant$Run == current_run), "Log2Intensity"] <- quant[which(quant$Run == current_run), "Log2Intensity"] + current_diff

}

# Recalculate intensity from log2intensity
quant$Intensity <- 2 ^ quant$Log2Intensity



# Process data with no normalization
processed <-
  dataProcess(
    quant,
    normalization = FALSE,
    summaryMethod = "TMP",
    cutoffCensored = "minFeature",
    censoredInt = "0",
    MBimpute = TRUE,
    maxQuantileforCensored = 0.999
  )

# Get subquant data
subquant <- quantification(processed)

# Write subquant to a file
write.table(x = subquant, file = paste(prefix, "-subquant.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

# Read in comparisons file
comparisons <- read.delim(file = "20230227-AK10-ph-comparisons.txt", quote = "")

# Grab condition columns from comparisons
comparisonmatrix <- comparisons[,2:ncol(comparisons)]

# Reorder conditions in alphabetical order (MSstast will fail if you don't do this)
comparisonmatrix <- as.matrix(comparisonmatrix[,order(names(comparisonmatrix))])

# Add names of comparisons as row names
row.names(comparisonmatrix) <- comparisons[,1]

# Run MSstats groupComparison function for comparisons in comparison file
compared <- groupComparison(contrast.matrix = comparisonmatrix, data = processed)

# Get results from MSstats comparison
results <- compared$ComparisonResult


write.table(x = results, file = paste(prefix, "-results.txt", sep = "", collapse = ""), sep = "\t", row.names = FALSE, quote = FALSE)


# Set log2fold-change and p-value cutoffs for downstream analyses
l2fc_max <- 1
l2fc_min <- -1
pval_max <- 0.05
adjpval_max <- 0.05

# Set database files
dbfiles <- ("20191010.SwissProt.Hsapiens.fasta")


# Get the number of databases
numdb <- length(dbfiles)

# Loop through each database file
for (i in 1:numdb) {
  # If it's the first database then read lines into dblines
  if (i == 1) {
    dblines <- readLines(con = dbfiles[i])
  } else {
    # If it's not the first database then read lines into tmp
    # Then concatenate tmp onto the end of dblines
    tmp <- readLines(con = dbfiles[i])
    dblines <- c(dblines, tmp)
  }
}

# Grab header lines starting with > character
headerlines <- dblines[grep(">", dblines)]

# Extract protein accession and protein name from header lines
protein.accession <- str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+.*)")[,2]
protein.description <- str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+.*)")[,3]
protein.name <- str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+)_.*")[,3]

# Create a named vector to lookup protein descriptions by accession
getProteinDescription <- protein.description
names(getProteinDescription) <- protein.accession
getProteinName <- protein.name
names(getProteinName) <- protein.accession

# Split accession column by semi-colons
split_acc <- str_split(string = results$Protein, pattern = ";")

# Go through row by row and extract protein accessions from site accessions
# Look up protein names and descriptions by lookup vectors
# Add Protein.Accession, Protein.Name, and Protein.Description columns to results table
# This section is very slow but I can't figure out a better way to do this
for (i in 1:length(split_acc)) {
  protein_accession <- str_match(split_acc[[i]], "([A-Z0-9]+)_ph.*")[,2]
  proteinAccessions <- paste(protein_accession, collapse = ";")
  results[i,"Protein.Accession"] <- proteinAccessions
  lookup_function <- function(x) unname(getProteinName[x])
  proteinNames <- paste(lapply(protein_accession, lookup_function), collapse = ";")
  results[i,"Protein.Name"] = proteinNames
  lookup_function <- function(x) unname(getProteinDescription[x])
  proteinDescriptions <- paste(lapply(protein_accession, lookup_function), collapse = ";")
  results[i,"Protein.Description"] <- proteinDescriptions
}

# Write results file to text
write.table(x = results, file = paste(prefix, "-results-ann.txt", sep = '',  collapse = ''), sep = "\t", quote = FALSE, row.names = FALSE)

# Melt and cast results to wide format and write to file
results_melted <- melt(data = results, id.vars = c("Protein", "Protein.Description", "Protein.Name", "Label"), measure.vars = c("log2FC", "pvalue", "adj.pvalue"))
results_wide <- dcast(data = results_melted, formula = Protein + Protein.Description + Protein.Name ~ Label + variable, value.var = "value")
write.table(x = results_wide, file = paste(prefix, "-results-wide.txt", collapse = '', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE)


# Write subquant data to file
write.table(x = subquant, file = paste(prefix, "-subquant.txt", sep = "", collapse = ""), sep = "\t", row.names = FALSE, quote = FALSE)



# Extract numerical data from subquant matrix for sample correlation heatmap
subquantmatrix <- as.matrix(subquant[,2:ncol(subquant)])

# Calculate pairwise correlations for all samples
cormat <- round(cor(subquantmatrix, use="complete.obs"),2)

# Row names are the same as column names
row.names(cormat) <- colnames(cormat)

# Set color scale from white to blue from 0 to 1 in 0.1 increments
breaks <- seq(0, 1, 0.1)
colorpalette <- colorRampPalette(c("white", "blue"))

# Plot heatmap of correlations and save to pdf
pdf(paste(prefix, "-sample-correlation-heatmap.pdf", collapse = '', sep = ''), width = 8, height = 8)
heatmap.2(x=cormat, col = colorpalette(100), trace="none", density.info = "none", margins = c(10,10), key.xlab = "Correlation")
dev.off()



# Extract numerical data from subquant matrix for sample correlation heatmap
subquantmatrix <- as.matrix(subquant[,2:ncol(subquant)])

# Add protein accession as row names
row.names(subquantmatrix) <- subquant$Protein

# Remove any rows with NA values (you will lose a lot of data!)
subquantmatrix <- subquantmatrix[which(!(is.na(rowSums(subquantmatrix)))),]

# Calculate principle components based on subquant matrix
pca <- prcomp(na.omit(subquantmatrix), center = TRUE, scale = TRUE)

# Combine pca x values and subquant data in a data frame
pca_subquant <- as.data.frame(cbind(pca$x, subquantmatrix))

# Add accession column
pca_subquant$Protein <- row.names(subquantmatrix)

# Write to a file
write.table(x = pca_subquant, file = "20230404-AK10-ph-subquant-pca.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Get the summary
pca_summary <- summary(pca)$importance

# Get the proportion of variance of PCs
pc1_propvar <- pca_summary[2,1]
pc2_propvar <- pca_summary[2,2]


# Extract PC1 and PC2 for plotting
pcaplot <- as.data.frame(cbind(pca$rotation[,1], pca$rotation[,2]))


# Parse drug from column names
timepoint <- str_match(string = colnames(subquantmatrix), pattern = "x*0*([^_]+)")[,2]


# Plot pca and color by cell type
ggplot(data=pcaplot, aes(x=pcaplot[,1], y=pcaplot[,2])) + geom_point(aes(color = timepoint), size = 3) + xlab(label = paste("Principle Component 1 ", pc1_propvar*100, sep = "")) + ylab(label = paste("Principle Component 2 ", pc2_propvar*100, sep = "")) + theme_bw() + theme(aspect.ratio = 1) + scale_color_discrete(name = "Time Point") + theme(aspect.ratio = 1)
ggsave(filename = paste(prefix, "-pca.pdf", sep = ""), width = 4, height = 4)





# Convert subquant to T/F based on NA values
subquant_tf <- is.na(subquant)[,2:ncol(subquant)]

# Sum each column
nacounts <- colSums(subquant_tf)

# Subtract na count from total count to get protien counts
proteincounts <- nrow(subquant_tf) - nacounts

# Convert to a data frame
proteincounts <- as.data.frame(proteincounts)

# Rename column of data frame
colnames(proteincounts) <- "Count"


# Add column of data frame to label Sample names
proteincounts$Sample <- row.names(proteincounts)

# Plot protein counts by samples
ggplot(data = proteincounts, aes(x = proteincounts$Sample, y = proteincounts$Count)) + geom_bar(stat = "identity", fill = "dodgerblue", color = "black") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(limits = c(0,15000), expand = c(0,0)) + xlab(label = "Sample") + ylab(label = "Num. Phosphosites Quantified")
ggsave(filename = paste(prefix, "-sitecounts-bysample.pdf", sep = ""), width = 6, height = 4)




#Melt subquant data to plot intensity box plots by sample
subquant_melted <- melt(data = subquant, id.vars = "Protein", measure.vars = colnames(subquant[,2:ncol(subquant)]), variable.name = "Sample", value.name = "Log2Intensity")


# Create a function to calculate the median of subquant log2intensity for each sample
subquantmedianfunc <- function (x) round(median(subquant[,x], na.rm = TRUE), 2)

# Apply function across each column of data in subquant
subquantmedians <- data.frame(unlist(lapply(X = 2:ncol(subquant), FUN = subquantmedianfunc)))

# Rename column header
colnames(subquantmedians) <- "Median Log2Intensity"

subquantmedians$Sample <- colnames(subquant[,2:ncol(subquant)])

# Plot

ggplot(data = subquant_melted, aes(x = subquant_melted$Sample, y = subquant_melted$Log2Intensity)) + geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90))  + xlab(label = NULL) + ylab(label = "Log2Intensity") + geom_text(data = subquantmedians, aes(x = subquantmedians$Sample, y = subquantmedians$`Median Log2Intensity`, label = subquantmedians$`Median Log2Intensity`), vjust = -0.5, size = 2)
ggsave(filename = paste(prefix, "-sample-intensities.pdf", sep = ""), width = 6, height = 4)




# Replace Inf/-Inf with NA values
results[which(results$log2FC == Inf), "log2FC"] <- NA
results[which(results$log2FC == -Inf), "log2FC"] <- NA


# Cast to wide format on l2fc
l2fc_wide <- dcast(data = results, formula = Protein + Protein.Name ~ Label, value.var = "log2FC")

# Get numerical data for correlation clustering
l2fc_mat <- as.matrix(l2fc_wide[, 3:ncol(l2fc_wide)])


# Calculate pairwise correlations for all samples
cormat <- round(cor(l2fc_mat, use="complete.obs"),2)

# Row names are the same as column names
row.names(cormat) <- colnames(cormat)

# Set color scale from white to blue from 0 to 1 in 0.1 increments
breaks <- seq(-1, 0, 1, 0.1)
colorpalette <- colorRampPalette(c("red", "white", "blue"))

# Plot heatmap of correlations and save to pdf
pdf(paste(prefix, "-l2fc-correlation-heatmap.pdf", collapse = '', sep = ''), width = 6, height = 6)
heatmap.2(x=cormat, col = colorpalette(100), trace="none", density.info = "none", margins = c(15,15), key.xlab = "Correlation")
dev.off()



# Set log2fold-change and p-value cutoffs for downstream analyses
l2fc_max <- 1
l2fc_min <- -1
pval_max <- 0.05
adjpval_max <- 0.05




# Extract differentially increased/decreased data
# Must be real numbers (i.e., not Inf, -Inf, or NA)
upregulated <-
  results[which((results$log2FC > l2fc_max) &
                  (results$pvalue < pval_max) &
                  (results$adj.pvalue < adjpval_max) &
                  (is.finite(results$log2FC)) &
                  (is.finite(results$pvalue))
  ), ]

downregulated <-
  results[which((results$log2FC < l2fc_min) &
                  (results$pvalue < pval_max) &
                  (results$adj.pvalue < adjpval_max) &
                  (is.finite(results$log2FC)) &
                  (is.finite(results$pvalue))
  ), ]

# Count number of up/down regulated sites using table function
upregulatedcounts <- as.data.frame(table(upregulated$Label))
downregulatedcounts <- as.data.frame(table(downregulated$Label))


# Label change direction so up and down regulated can be combined in a single data frame
upregulatedcounts$Direction <- "Increased"
downregulatedcounts$Direction <- "Decreased"
significantcounts <- rbind(upregulatedcounts, downregulatedcounts)


# Flip the sign on decreased proteins
significantcounts[which(significantcounts$Direction == "Decreased"), "Freq"] <- -1 * significantcounts[which(significantcounts$Direction == "Decreased"), "Freq"]

# Make a column for vjust
significantcounts$vjust <- -0.5
significantcounts[which(significantcounts$Direction == "Decreased"), "vjust"] <- 1.5

# Get time points from Var1
significantcounts$TimePoint <- str_match(string = significantcounts$Var1, pattern = "PIC_(\\S+)-Control")[,2]

# Set colors for bar plots
colors <- c("Increased" = "blue", "Decreased" = "red")


ggplot(data = significantcounts, aes(x = significantcounts$TimePoint, y = significantcounts$Freq)) + geom_bar(stat = "identity", aes(fill = factor(significantcounts$Direction))) + scale_fill_manual(values = colors, name = "Direction") + geom_text(label = abs(significantcounts$Freq), vjust = significantcounts$vjust) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + xlab(label = "Time Point") + ylab(label = "Num. Sig. Phospho. Changes") + ggtitle(label = "|Log2FC| > 1 & adj.pval < 0.05") + scale_x_discrete(limits = c("05m", "30m", "02h", "24h")) + scale_y_continuous(limits = c(-2000, 2500)) + theme(aspect.ratio = 1)
ggsave(filename = paste(prefix, "-diff-abun-barplots-adjpval05.pdf", sep = ""), width = 4, height = 4)



# Read in abundance results
results_ab <- read.delim(file = "../../20230227-AK10-ab//msstats/20230227-AK10-results-ann.txt", quote = "", stringsAsFactors = FALSE)

# Make a lookup table for protein abundance and adj.pval
getAbLog2FC <- results_ab$log2FC
names(getAbLog2FC) <- paste(results_ab$Protein, results_ab$Label)
getAbAdjPval <- results_ab$adj.pvalue
names(getAbAdjPval) <- paste(results_ab$Protein, results_ab$Label)

# Add abundance log2fc to ph results
results$log2FC_abun <- getAbLog2FC[paste(results$Protein.Accession, results$Label)]
results$adj.pvalue_abun <- getAbAdjPval[paste(results$Protein.Accession, results$Label)]

# Get a list of comparisons
comparisons <- unique(results$Label)


# Loop through comparison
for (i in 1:length(comparisons)) {

  # Get current comparison
  current_comparison <- comparisons[i]

  # Filter results for current comparison
  current_results <- results[which(results$Label == current_comparison),]

  # Plot l2fc ab vs ph
  ggplot(data = current_results, aes(x = current_results$log2FC_abun, y = current_results$log2FC)) + geom_point(size = 0.5) + stat_density2d(bins = 30, color = "red") + theme_bw() + xlab(label = "Log2FC Protein Abundance") + ylab(label = "Log2FC Phosphorylation") + ggtitle(label = current_comparison) + scale_x_continuous(limits = c(-6,6)) + scale_y_continuous(limits = c(-6,6))
  ggsave(filename = paste("20230323-AK10-ph-v-ab-", current_comparison, ".pdf", sep = ""), width = 4, height = 4)

}



# Make volcano plots of ph with fill color by ab log2fc
# Loop through comparison
for (i in 1:length(comparisons)) {

  # Get current comparison
  current_comparison <- comparisons[i]

  # Filter results for current comparison
  current_results <- results[which(results$Label == current_comparison),]

  # Filter out NA values in l2fc or pvalue
  current_results <- current_results[which(!(is.na(current_results$log2FC))),]
  current_results <- current_results[which(!(is.na(current_results$pvalue))),]

  # Add a logp column
  current_results$logP <- -1 * log(current_results$pvalue) / log(10)

  ggplot(data = current_results, aes(x = current_results$log2FC, y = current_results$logP)) + geom_point(aes(color = current_results$log2FC_abun)) + scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "Ab Log2FC") + theme_bw() + xlab(label = "Log2FC Phosphorylation") + ylab(label = "P-value Phosphorylation") + ggtitle(label = current_comparison) + theme(aspect.ratio = 1)
  ggsave(filename = paste("20230323-AK10-ph-volcano-fillbyab-", current_comparison, ".pdf", sep = ""), width = 6, height = 4)

}


# Count sig changes in ph that also are sig changes in ab

# Make a table to store data
sig_ab_ph_table <- data.frame(
  comparison = character(),
  nsig_ph = integer(),          # Num sig ph sites all
  sig_ph_sig_ab = integer(),    # Num sig ph sites also sig ab
  sig_ph_quant_ab = integer()   # Num sig ph sites all quant ab
)

# Loop through comparisons
for (i in 1:length(comparisons)) {

  # Get Current comparison
  current_comparison <- comparisons[i]

  # Count ph sig changes regardless of ab data
  nsig_ph <- nrow(results[which(
    (results$Label == current_comparison) &
      (is.finite(results$log2FC)) &
      (abs(results$log2FC) > 1) &
      (results$adj.pvalue < 0.05)),])


  # Count sig changes ph with finite ab l2fc and pval
  sig_ph_quant_ab <- nrow(results[which(
    (results$Label == current_comparison) &
      (is.finite(results$log2FC)) &
      (abs(results$log2FC) > 1) &
      (results$adj.pvalue < 0.05) &
      (is.finite(results$log2FC_abun)) &
      (is.finite(results$adj.pvalue_abun))),])

  # Count sig changes ph and ab
  sig_ph_sig_ab <- nrow(results[which(
    (results$Label == current_comparison) &
      (is.finite(results$log2FC)) &
      (abs(results$log2FC) > 1) &
      (results$adj.pvalue < 0.05) &
      (is.finite(results$log2FC_abun)) &
      (abs(results$log2FC_abun > 1)) &
      (results$adj.pvalue_abun < 0.05)),])

  sig_ab_ph_table <- rbind(data.frame(
    comparison = current_comparison,
    nsig_ph = nsig_ph,
    sig_ph_sig_ab = sig_ph_sig_ab,
    sig_ph_quant_ab = sig_ph_quant_ab), sig_ab_ph_table)


}

# Write table to file
write.table(x = sig_ab_ph_table, file = "20230323-AK10-ph-sigph-sigab.txt", sep = "\t", quote = FALSE, row.names = FALSE)






# Cast l2fc to wide format from results
l2fc_wide <- dcast(data = results, formula = Protein + Protein.Name ~ Label, value.var = "log2FC")

# Get proteins with non-finite row sums
l2fc_wide <- l2fc_wide[which(is.finite(rowSums(l2fc_wide[,3:ncol(l2fc_wide)]))),]

# Get numeric data as a matrix
l2fc_mat <- l2fc_wide[,3:ncol(l2fc_wide)]

# Calculate principle components based on subquant matrix
pca <- prcomp(na.omit(l2fc_mat), center = TRUE, scale = TRUE)


# Get importance table from summary
importance <- summary(pca)$importance

# Get proportion of variances for each component
propvar <- as.data.frame(importance[2,])

# Rename column
colnames(propvar) <- "Proportion_variance"


# Extract PC1 and PC2 for plotting
pcaplot <- as.data.frame(cbind(pca$rotation[,1], pca$rotation[,2])

# Create a column for time point
pcaplot$TimePoint <- str_match(string = row.names(pcaplot), pattern = "PIC_(\\S+)-Control")[,2]


# Plot pc1 vs pc2
ggplot(data = pcaplot, aes(x = pcaplot$V1, y = pcaplot$V2)) + geom_point(aes(color = pcaplot$TimePoint), size = 3) + theme_bw() + xlab(label = paste("Principal Component 1", propvar["PC1",]*100)) + ylab(label = paste("Principal Component 2", propvar["PC2",]*100)) + scale_color_discrete(name = "Time Point") + theme(aspect.ratio = 1)
ggsave(filename = "20230404-AK10-ph-l2fc-pca.pdf", width = 4, height = 4)
