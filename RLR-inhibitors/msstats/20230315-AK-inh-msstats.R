library(MSstats)
library(reshape2)
library(dplyr)
library(ggplot2)
library(gplots)
library(stringr)

# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital/Data/JohnsonLab/AndrewKurland/PolyIC-Inhibitor-Screen/20230315-AKinh-analysis/20230315-AKinh-msstats/")



# Set prefix for naming output files
prefix <- "20230405-AK-inh-ab"

# Set log2fold-change and p-value cutoffs for downstream analyses
l2fc_max <- 1
l2fc_min <- -1
pval_max <- 0.05
adjpval_max <- 0.05


# Set database files
dbfiles <- ("20191010.SwissProt.Hsapiens.fasta")


# Read in Spectronaut report
raw <- read.delim(file = "20221025_124220_20221021-AKinh_Report.tsv", quote = "")


# Drop 1 uM conditions
raw <- raw[which(!(str_detect(string = raw$R.Condition, pattern = "_1uM"))),]

# Convert to MSstats format
quant <- SpectronauttoMSstatsFormat(input = raw, intensity = "NormalizedPeakArea")


# Run MSstats dataProcess function with no normalization
processed <- dataProcess(quant, normalization = FALSE, summaryMethod = "TMP", censoredInt = "0", MBimpute = TRUE, maxQuantileforCensored = 0.999, featureSubset = "topN", n_top_feature = 10)


# Get subquant data from processed
subquant <- quantification(processed)

# Write subquant file to text
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
pdf(paste(prefix, "-sample-correlation-heatmap.pdf", collapse = '', sep = ''), width = 16, height = 16)
heatmap.2(x=cormat, col = colorpalette(100), trace="none", density.info = "none", margins = c(10,10), key.xlab = "Correlation")
dev.off()


# Remove MRT50 columns
subquant <- subquant[,which(!str_detect(string = colnames(subquant), pattern = "MRT50"))]

# Extract numerical data from subquant matrix for sample correlation heatmap
subquantmatrix <- as.matrix(subquant[,2:ncol(subquant)])



# Calculate principle components based on subquant matrix
pca <- prcomp(na.omit(subquantmatrix), center = TRUE, scale = TRUE)

# Extract PC1 and PC2 for plotting
pcaplot <- as.data.frame(cbind(pca$rotation[,1], pca$rotation[,2]))

# Get condition from sample name
pcaplot$Condition <- str_match(string = row.names(pcaplot), pattern = "(\\S+)\\s*_\\d+")[,2]

# Create a column for conditions to label with text and fill with NAs
pcaplot$LabeledCondition <- NA

# Label points for hits that separate
pcaplot[which(pcaplot$Condition == "MRT10"), "LabeledCondition"] <- "TBK1 low"
pcaplot[which(pcaplot$Condition == "MRT20."), "LabeledCondition"] <- "TBK1 high"
pcaplot[which(pcaplot$Condition == "AMPKA1_10uM"), "LabeledCondition"] <- "AMPKA1"
pcaplot[which(pcaplot$Condition == "PLK1_10uM"), "LabeledCondition"] <- "PLK1"
pcaplot[which(pcaplot$Condition == "ATR_10uM"), "LabeledCondition"] <- "ATR"
pcaplot[which(pcaplot$Condition == "DMSO"), "LabeledCondition"] <- "DMSO"



# Plot pca and color by cell type
ggplot(data=pcaplot, aes(x=pcaplot[,1], y=pcaplot[,2])) + geom_point(aes(color = pcaplot$LabeledCondition), size = 3) + xlab(label = "Principle Component 1") + ylab(label = "Principle Component 2") + theme_bw() + theme(aspect.ratio = 1) + scale_color_discrete(name = "Condition") + geom_text(aes(label = pcaplot$LabeledCondition, color = pcaplot$LabeledCondition), hjust = -0.1, vjust = 0.1)
ggsave(filename = paste(prefix, "-pca-colors.pdf", sep = ""), width = 5, height = 5)


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
ggplot(data = proteincounts, aes(x = proteincounts$Sample, y = proteincounts$Count)) + geom_bar(stat = "identity", fill = "dodgerblue", color = "black") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(limits = c(0,6000), expand = c(0,0)) + xlab(label = "Sample") + ylab(label = "Num. Proteins Quantified")
ggsave(filename = paste(prefix, "-proteincounts-bysample.pdf", sep = ""), width = 6, height = 4)


#Melt subquant data to plot intensity box plots by sample
subquant_melted <- melt(data = subquant, id.vars = "Protein", measure.vars = colnames(subquant[,2:ncol(subquant)]), variable.name = "Sample", value.name = "Log2Intensity")


# Aggregate by medians
subquantmedians <- aggregate(Log2Intensity ~ Sample, subquant_melted, median)


# Plot

ggplot(data = subquant_melted, aes(x = subquant_melted$Sample, y = subquant_melted$Log2Intensity)) + geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90))  + xlab(label = NULL) + ylab(label = "Log2Intensity") + geom_text(data = subquantmedians, aes(x = subquantmedians$Sample, y = subquantmedians$Log2Intensity, label = round(subquantmedians$Log2Intensity, digits = 1)), vjust = -0.5, size = 2)
ggsave(filename = paste(prefix, "-sample-intensities.pdf", sep = ""), width = 16, height = 6)



# Read in comparisons file
comparisons <- read.delim(file = "Comparison_AKinh.txt", quote = "")

# Grab condition columns from comparisons
comparisonmatrix <- comparisons[,2:ncol(comparisons)]

# Reorder conditions in alphabetical order (MSstast will fail if you don't do this)
comparisonmatrix <- as.matrix(comparisonmatrix[,order(names(comparisonmatrix))])

# Add names of comparisons as row names
row.names(comparisonmatrix) <- comparisons[,1]

# Reorder columns for MSstats (capital letters first)
comparisonmatrix <- comparisonmatrix[, c(1,2,3,4,5,6,7,8,9,10,11,13,14,15,12,19,16,17,18)]

# MRT20 has two spaces after it >:(
colnames(comparisonmatrix)[13] <- "MRT20 "

# Run comparison function in MSstats
compared <- groupComparison(contrast.matrix = comparisonmatrix, data = processed)


# Store MSstats results (log2fold-changes and adjusted p-values) in results data frame
results <- compared$ComparisonResult

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
# This section is a little slow but I can't figure out a better way to do this
for (i in 1:length(split_acc)) {
  protein_accession <- str_match(split_acc[[i]], "([A-Z0-9]+).*")[,2]
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

# Drop MRT50 from resutls (cells are not viable) and mock_DMSO
results <- results[which(!(results$Label == "MRT50_DMSO")),]
results <- results[which(!(results$Label == "Mock_DMSO")),]

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


# Set colors for bar plots
colors <- c("Increased" = "red", "Decreased" = "blue")

# Remove MRT50
significantcounts <- significantcounts[which(!(significantcounts$Var1 == "MRT50_DMSO")),]

# Replace labels for plotting
significantcounts[which(significantcounts$Var1 == "AMPKA1_10uM_DMSO"), "PlotLabel"] <- "AMPKA1 / DMSO"
significantcounts[which(significantcounts$Var1 == "ATM_10uM_DMSO"), "PlotLabel"] <- "ATM / DMSO"
significantcounts[which(significantcounts$Var1 == "ATR_10uM_DMSO"), "PlotLabel"] <- "ATR / DMSO"
significantcounts[which(significantcounts$Var1 == "AURKA_10uM_DMSO"), "PlotLabel"] <- "AURKA / DMSO"
significantcounts[which(significantcounts$Var1 == "AURKB_10uM_DMSO"), "PlotLabel"] <- "AURKB / DMSO"
significantcounts[which(significantcounts$Var1 == "DMSO_Control"), "PlotLabel"] <- "Z DMSO / Neg. Control"
significantcounts[which(significantcounts$Var1 == "ERK12_10uM_DMSO"), "PlotLabel"] <- "ERK12 / DMSO"
significantcounts[which(significantcounts$Var1 == "ERK2_10uM_DMSO"), "PlotLabel"] <- "ERK2 / DMSO"
significantcounts[which(significantcounts$Var1 == "GSK3AB_10uM_DMSO"), "PlotLabel"] <- "GSK3AB / DMSO"
significantcounts[which(significantcounts$Var1 == "JNK123_10uM_DMSO"), "PlotLabel"] <- "JNK123 / DMSO"
significantcounts[which(significantcounts$Var1 == "MRT10_DMSO"), "PlotLabel"] <- "TBK1 Low / DMSO"
significantcounts[which(significantcounts$Var1 == "MRT20 _DMSO"), "PlotLabel"] <- "TBK1 High / DMSO"
significantcounts[which(significantcounts$Var1 == "mTOR_10uM_DMSO"), "PlotLabel"] <- "mTOR / DMSO"
significantcounts[which(significantcounts$Var1 == "p38a_10uM_DMSO"), "PlotLabel"] <- "p38a / DMSO"
significantcounts[which(significantcounts$Var1 == "p38ab_10uM_DMSO"), "PlotLabel"] <- "p38ab / DMSO"
significantcounts[which(significantcounts$Var1 == "PLK1_10uM_DMSO"), "PlotLabel"] <- "PLK1 / DMSO"




ggplot(data = significantcounts, aes(x = significantcounts$PlotLabel, y = significantcounts$Freq)) + geom_bar(stat = "identity", aes(fill = factor(significantcounts$Direction))) + scale_fill_manual(values = colors, name = "Direction") + geom_text(label = abs(significantcounts$Freq), vjust = significantcounts$vjust) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + xlab(label = "Comparison") + ylab(label = "Num. Sig. Changes") + ggtitle(label = "|Log2FC| > 1 & adj.pval < 0.05") + scale_y_continuous(limits = c(-500,200))
ggsave(filename = paste(prefix, "-diff-abun-barplots-adjpval05.pdf", sep = ""), width = 6, height = 4)




# Replace labels in results for plotting
results[which(results$Label == "AMPKA1_10uM_DMSO"), "PlotLabel"] <- "AMPKA1 / DMSO"
results[which(results$Label == "ATM_10uM_DMSO"), "PlotLabel"] <- "ATM / DMSO"
results[which(results$Label == "ATR_10uM_DMSO"), "PlotLabel"] <- "ATR / DMSO"
results[which(results$Label == "AURKA_10uM_DMSO"), "PlotLabel"] <- "AURKA / DMSO"
results[which(results$Label == "AURKB_10uM_DMSO"), "PlotLabel"] <- "AURKB / DMSO"
results[which(results$Label == "DMSO_Control"), "PlotLabel"] <- "Z DMSO / Neg. Control"
results[which(results$Label == "ERK12_10uM_DMSO"), "PlotLabel"] <- "ERK12 / DMSO"
results[which(results$Label == "ERK2_10uM_DMSO"), "PlotLabel"] <- "ERK2 / DMSO"
results[which(results$Label == "GSK3AB_10uM_DMSO"), "PlotLabel"] <- "GSK3AB / DMSO"
results[which(results$Label == "JNK123_10uM_DMSO"), "PlotLabel"] <- "JNK123 / DMSO"
results[which(results$Label == "MRT10_DMSO"), "PlotLabel"] <- "TBK1 Low / DMSO"
results[which(results$Label == "MRT20 _DMSO"), "PlotLabel"] <- "TBK1 High / DMSO"
results[which(results$Label == "mTOR_10uM_DMSO"), "PlotLabel"] <- "mTOR / DMSO"
results[which(results$Label == "p38a_10uM_DMSO"), "PlotLabel"] <- "p38a / DMSO"
results[which(results$Label == "p38ab_10uM_DMSO"), "PlotLabel"] <- "p38ab / DMSO"
results[which(results$Label == "PLK1_10uM_DMSO"), "PlotLabel"] <- "PLK1 / DMSO"

# Read in a file of ISGs
isgs <- read.delim(file = "ISGs.txt", quote = "", stringsAsFactors = FALSE)

# Make a lookup table for isg gene names
getGeneName <- isgs$Gene.Names..primary.
names(getGeneName) <- isgs$Entry

# Filter results for ISGs
results_isgs <- results[which(results$Protein %in% isgs$Entry),]

# Add gene names to results
results_isgs$Gene.Name <- getGeneName[as.character(results_isgs$Protein)]

# Get isgs with l2fc > 1 and adj.pval < 0.05 in DMSO vs Control
sig_isgs <- unique(results_isgs[which(
  (results_isgs$log2FC > 1) &
    (results_isgs$adj.pvalue < 0.05) &
    (results_isgs$Label %in% "DMSO_Control")), "Protein"])


# Filter for sig isgs
results_isgs <- results_isgs[which(results_isgs$Protein %in% sig_isgs),]


# Cast to wide format on log2fc
results_isgs_wide <- dcast(data = results_isgs, formula = Protein + Protein.Name + Gene.Name ~ Label, value.var = "log2FC")

# FIlter out rows with non-finite rowSums
results_isgs_wide <- results_isgs_wide[which(is.finite(rowSums(results_isgs_wide[,4:ncol(results_isgs_wide)]))),]

# Get numeric data for heatmap plotting
l2fc_mat <- as.matrix(results_isgs_wide[, 4:ncol(results_isgs_wide)])

# Make row names protein names
row.names(l2fc_mat) <- results_isgs_wide$Gene.Name

# Set color scale from red to white to blue from 0 to 1 in 0.1 increments
breaks <- seq(-1, 1, 0.1)
colorpalette <- colorRampPalette(c("blue", "white", "red"))

# Plot heatmap of correlations and save to pdf
pdf(paste(prefix, "-ISG-L2FC-heatmap.pdf", collapse = '', sep = ''), width = 8, height = 5)
heatmap.2(x=l2fc_mat, col = colorpalette(100), trace="none", density.info = "none", margins = c(10,10), key.xlab = "Log2FC")
dev.off()





# Filter subquant for ISGs
subquant_isgs <- subquant[which(subquant$Protein %in% isgs$Entry),]

# Melt down
subquant_isgs_melted <- melt(data = subquant_isgs, id.vars = "Protein", variable.name = "Sample", value.name = "Log2Intensity")

# Add gene names
subquant_isgs_melted$Gene.Name <- getGeneName[as.character(subquant_isgs_melted$Protein)]

# Parse out condition from sample name
subquant_isgs_melted$Condition <- str_match(string = subquant_isgs_melted$Sample, pattern = "(\\S+)\\s*_\\d+")[,2]

# Aggregate by medians
subquant_isgs_melted_agg <- aggregate(Log2Intensity ~ Protein + Gene.Name + Condition, subquant_isgs_melted, median)

# Cast to wide format
subquant_isgs_wide <- dcast(data = subquant_isgs_melted_agg, formula = Protein + Gene.Name ~ Condition, value.var = "Log2Intensity")

# Get numeric data for heatmap plotting
l2int_mat <- as.matrix(subquant_isgs_wide[, 3:ncol(subquant_isgs_wide)])

# Replace NAs with zeros
l2int_mat[is.na(l2int_mat)] <- 0

# Name rows with gene name
row.names(l2int_mat) <- subquant_isgs_wide$Gene.Name

# Set color scale from blue to white to red from 0 to 1 in 0.1 increments
breaks <- seq(0, 20, 0.1)
colorpalette <- colorRampPalette(c("white", "white", "white", "red"))

# Plot heatmap of correlations and save to pdf
pdf(paste(prefix, "-ISG-Log2Intensity-heatmap.pdf", collapse = '', sep = ''), width = 10, height = 8)
heatmap.2(x=l2int_mat, col = colorpalette(100), trace="none", density.info = "none", margins = c(10,10), key.xlab = "Log2Intensity")
dev.off()


# Remove Mock_DMSO comparison and MRT50-Control
results <- results[which(!(results$Label == "Mock_DMSO")),]
results <- results[which(!(results$Label == "MRT50_Control")),]

# Fix MRT20 label
results[which(results$Label == "MRT20 _Control"), "Label"] <- "MRT20_Control"

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
pcaplot <- as.data.frame(cbind(pca$rotation[,1], pca$rotation[,2]))



# Get comparison name from row names
pcaplot$Comparison <- row.names(pcaplot)

# Get drug target from comparison
pcaplot$Target <- str_match(string = pcaplot$Comparison, pattern = "(\\S+)_Control")[,2]

# Remove 10uM text from target
pcaplot$Target <- str_replace(string = pcaplot$Target, pattern = "_10uM", replacement = "")

# Plot pc1 vs pc2
ggplot(data = pcaplot, aes(x = pcaplot$V1, y = pcaplot$V2)) + geom_point(size = 3) + theme_bw() + xlab(label = paste("Principal Component 1", propvar["PC1",]*100)) + ylab(label = paste("Principal Component 2", propvar["PC2",]*100)) + scale_color_discrete(name = "Time Point") + theme(aspect.ratio = 1) + geom_text(aes(x = pcaplot$V1, y = pcaplot$V2, label = pcaplot$Target), hjust = -0.25)
ggsave(filename = "20230405-AKinh-l2fc-pca.pdf", width = 4, height = 4)



