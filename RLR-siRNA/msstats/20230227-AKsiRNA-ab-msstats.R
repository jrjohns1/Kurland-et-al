library(MSstats)
library(reshape2)
library(dplyr)
library(ggplot2)
library(gplots)
library(stringr)


# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital/Data/JohnsonLab/AndrewKurland/PolyIC-siRNA/msstats/")


# Set prefix for naming output files
prefix <- "20230327-AKsiRNA"

# Set log2fold-change and p-value cutoffs for downstream analyses
l2fc_max <- 1
l2fc_min <- -1
pval_max <- 0.05
adjpval_max <- 0.05

# Set database files
dbfiles <- ("~/OneDrive - The Mount Sinai Hospital/DB/20191010.SwissProt.Hsapiens.fasta")

# Read in Spectronaut report
raw <- read.delim(file = "20230327_124255_0230327-AksiRNA_Report.tsv", quote = "")

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
pdf(paste(prefix, "-sample-correlation-heatmap.pdf", collapse = '', sep = ''), width = 8, height = 8)
heatmap.2(x=cormat, col = colorpalette(100), trace="none", density.info = "none", margins = c(10,10))
dev.off()


# Extract numerical data from subquant matrix for sample correlation heatmap
subquantmatrix <- as.matrix(subquant[,2:ncol(subquant)])

# Calculate principle components based on subquant matrix
pca <- prcomp(na.omit(subquantmatrix), center = TRUE, scale = TRUE)

# Extract PC1 and PC2 for plotting
pcaplot <- as.data.frame(cbind(pca$rotation[,1], pca$rotation[,2]))

# Get condition from row names
pcaplot$condition <- str_match(string = row.names(pcaplot), pattern = "(\\S+)_\\d")[,2]

# Plot pca and color by cell type
ggplot(data=pcaplot, aes(x=pcaplot[,1], y=pcaplot[,2])) + geom_point(aes(color = pcaplot$condition), size = 3) + xlab(label = "Principle Component 1") + ylab(label = "Principle Component 2") + theme_bw() + theme(aspect.ratio = 1) + scale_color_discrete(name = "Time Point") + theme(aspect.ratio = 1)
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
ggplot(data = proteincounts, aes(x = proteincounts$Sample, y = proteincounts$Count)) + geom_bar(stat = "identity", fill = "dodgerblue", color = "black") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(limits = c(0,4750), expand = c(0,0)) + xlab(label = "Sample") + ylab(label = "Num. Proteins Quantified")
ggsave(filename = paste(prefix, "-proteincounts-bysample.pdf", sep = ""), width = 6, height = 4)


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







# Read in comparisons file
comparisons <- read.delim(file = "20230327-AKsiRNA-comparisons.txt", quote = "")

# Grab condition columns from comparisons
comparisonmatrix <- comparisons[,2:ncol(comparisons)]

# Reorder conditions in alphabetical order (MSstast will fail if you don't do this)
comparisonmatrix <- as.matrix(comparisonmatrix[,order(names(comparisonmatrix))])

# Add names of comparisons as row names
row.names(comparisonmatrix) <- comparisons[,1]

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
breaks <- seq(1, 1, 0.1)
colorpalette <- colorRampPalette(c("blue", "white", "red"))

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


# Set colors for bar plots
colors <- c("Increased" = "red", "Decreased" = "blue")


ggplot(data = significantcounts, aes(x = significantcounts$Var1, y = significantcounts$Freq)) + geom_bar(stat = "identity", aes(fill = factor(significantcounts$Direction))) + scale_fill_manual(values = colors, name = "Direction") + geom_text(label = abs(significantcounts$Freq), vjust = significantcounts$vjust) + theme_bw()  + xlab(label = "Comparison") + ylab(label = "Num. Sig. Protein Changes") + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(limits = c(-300,300))
ggsave(filename = paste(prefix, "-diff-abun-barplots-adjpval05.pdf", sep = ""), width = 6, height = 4)



# Get a list of ISGs
isgs <- read.delim(file = "~/OneDrive - The Mount Sinai Hospital/DB/ISGs.txt", quote = "", stringsAsFactors = FALSE)





# Filter results for isgs
results_isgs <- results[which(results$Protein %in% isgs$Entry),]

# Get isgs with l2fc > 1 and adj.pval < 0.05 in siControl_INF or siControl_PI vs siControl
sig_isgs <- unique(results_isgs[which(
  (results_isgs$log2FC > 1) &
    (results_isgs$adj.pvalue < 0.05) &
    (results_isgs$Label %in% c("siControl_INF-siControl", "siControl_PIC-siControl"))), "Protein"])


# FIlter for sig isgs
results_isgs <- results_isgs[which(results_isgs$Protein %in% sig_isgs),]

# Parse out whether comparisons are PIC, IFN, or no treatment
results_isgs$Treatment <- "None"
results_isgs[which(str_detect(string = results_isgs$Label, pattern = "_PIC")), "Treatment"] <- "PIC"
results_isgs[which(str_detect(string = results_isgs$Label, pattern = "_INF")), "Treatment"] <- "IFN"

# Cast to wide format on log2fc
results_isgs_wide <- dcast(data = results_isgs, formula = Protein + Protein.Name ~ Label, value.var = "log2FC")

# Get numeric data for clustering
l2fc_mat <- as.matrix(results_isgs_wide[,3:ncol(results_isgs_wide)])

# Replace NAs with 0s
l2fc_mat[is.na(l2fc_mat)] <- 0


# Calculate distances
d <- dist(l2fc_mat)

# Cluster
clustered <- hclust(d)


# Plot as clustered heatmap
ggplot(data = results_isgs, aes(x = results_isgs$Label, y = results_isgs$Protein.Name)) + geom_tile(aes(fill = results_isgs$log2FC)) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Log2FC") + scale_y_discrete(limits = results_isgs_wide[clustered$order, "Protein.Name"], expand = c(0,0)) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(expand = c(0,0)) + xlab(label = "Comparison") + ylab(label = "Protein") + facet_wrap(Treatment ~ ., scales = "free")
ggsave(filename = paste(prefix, "-ISG-heatmaps.pdf", sep = ""), width = 6, height = 6)




# Plot siATM_PIC vs siControl_PIC as a scatterplot
ggplot(data = results_isgs_wide, aes(x = results_isgs_wide$`siATM_PIC-siControl`, y = results_isgs_wide$`siControl_PIC-siControl`)) + geom_point()

t.test(x = results_isgs_wide$`siAMPKA1_PIC-siControl`, y = results_isgs_wide$`siControl_PIC-siControl`)
