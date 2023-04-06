library(ksea)
library(reshape2)
library(ggplot2)
library(stringr)
library(readxl)


# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital/Data/JohnsonLab/AndrewKurland/PolyIC-Timecourse/20230227-AK10-analysis/20230227-AK10-ph/ksea/")



# Load ProtMapper kinase substrate table
protmapper <- as.data.frame(read_excel(path = "~/OneDrive - The Mount Sinai Hospital/DB/20230322-ProtMapper-kinase-substrate-table.xlsx"))

# Construct a phacc column
protmapper$PH_ACC <- paste(protmapper$TARGET_UP_ID, "_ph", protmapper$TARGET_POS, sep = "")

# Filter for belief = 1
protmapper <- protmapper[which(protmapper$BELIEF == 1),]

# Filter for kinases
protmapper <- protmapper[which(protmapper$CTRL_IS_KINASE == TRUE),]

# Make a table of kinases and substrates
kinsub <- protmapper[, c("CTRL_ID", "PH_ACC")]

# Change column names
colnames(kinsub) <- c("KINASE", "SUB_PH_ACC")

# Create a list to store kinase-substrate annotations
kinases_gsea <- list()

# Get list of kinases
kinases <- unique(kinsub$KINASE)

# Make a lookup table from kinase name to gene name
getGeneName <- protmapper$CTRL_GENE_NAME
names(getGeneName) <- protmapper$CTRL_ID


# Loop through each kinase and load into the list
for (i in 1:length(kinases)) {

  # Get the substrates of the current kinase
  #substrates <- unique(kinsub_allph[which(kinsub_allph$KINASE == kinases[i]), "PhAcc"])

  substrates <- unique(kinsub[which(kinsub$KINASE == kinases[i]), "SUB_PH_ACC"])

  # Store list of substrates in the kinases_gsea list
  kinases_gsea[[i]] <- substrates

}

# Name each element in kinases_gsea list with KINASE name
names(kinases_gsea) <- kinases






# Read in results
results <- read.delim(file = "../msstats/20230227-AK10-ph-results-ann.txt", quote = "", stringsAsFactors = FALSE)

# Get a list of comparisons
comparisons <- unique(results$Label)

# Create a data frame to store GSEA results
all_ksea_results <- data.frame(
  comparison = character(),
  kinase = character(),
  ES = double(),
  pval = double(),
  nsub = double()
)


# Loop through comparisons
for (i in 1:length(comparisons)) {

  # Get current comparison
  current_comparison <- comparisons[i]

  # Print current comparison
  print(current_comparison)

  # Get l2fc data for current comparison
  current_l2fc <- results[which((results$Label == current_comparison) & (is.finite(results$log2FC))), "log2FC"]

  # Get names
  current_l2fc_names <- results[which((results$Label == current_comparison) & (is.finite(results$log2FC))), "Protein"]

  # Rank order l2fc
  current_l2fc_ranked <- current_l2fc[order(current_l2fc, decreasing = TRUE)]

  # Order corresponding protein names by same ranking
  names(current_l2fc_ranked) <- current_l2fc_names[order(current_l2fc, decreasing = TRUE)]


  #ksea_result <- ksea_batchKinases(ranking = names(current_l2fc_ranked), logvalues = current_l2fc_ranked, regulons = kinases_gsea, trial = 1000)

  # Loop through each kinase
  for (j in 1:length(kinases)) {

    # Get the current kinase
    current_kinase <- kinases[j]

    # Print the current kinase
    print(current_kinase)

    # Check that we have measurements for at least 3 substrates

    nsub <- length(intersect(names(current_l2fc_ranked), kinases_gsea[[current_kinase]]))

    # Run GSEA for this comparison if nsub > 3
    if (nsub > 3) {

      # Run ksea
      ksea_result <- ksea(names(current_l2fc_ranked), current_l2fc_ranked, kinases_gsea[[current_kinase]], trial=1000, significance = TRUE)

      # Print ksea result
      print (paste(current_comparison, j, "of", length(kinases), current_kinase, ksea_result$ES, ksea_result$p.value))

      # Add ksea result to table of all results
      all_ksea_results <- rbind(
        data.frame(
          comparison = current_comparison,
          kinase = current_kinase,
          ES = ksea_result$ES,
          pval = ksea_result$p.value,
          nsub = nsub
        ), all_ksea_results
      )

    }
  }
}


# Add gene name to ksea results
all_ksea_results$kinase_gene <- getGeneName[all_ksea_results$kinase]

# Adjust p-values
all_ksea_results$adj.pval <- p.adjust(p = all_ksea_results$pval, method = "fdr")

# Write to a file
write.table(x = all_ksea_results, file = "20230322-AK10-ph-ksea.txt", sep = "\t", row.names = FALSE, quote = FALSE)




# Get kinases with adj.pval < 0.05 in any comparison
sig_kinases <- unique(all_ksea_results[which(all_ksea_results$adj.pval < 0.05), "kinase"])

# Filter ksea results for sig kinases
sig_ksea <- all_ksea_results[which(all_ksea_results$kinase %in% sig_kinases),]


# Add logp column
sig_ksea$logp <- -1 * log(sig_ksea$pval) / log(10)

# Convert inf logp to 3
sig_ksea[which(sig_ksea$logp == Inf), "logp"] <- 3

# Remove kinases with NA (these are all duplicated)
sig_ksea <- sig_ksea[which(!(is.na(sig_ksea$kinase_gene))),]

# Remove P17252 because it is bovine origin
sig_ksea <- sig_ksea[which(!(sig_ksea$kinase == "P17252")),]

# Flip the sign on logp for negative ES values
sig_ksea[which(sig_ksea$ES < 0), "logp"] <- -1 * sig_ksea[which(sig_ksea$ES < 0), "logp"]

# Make all ES values positive
sig_ksea$ES <- abs(sig_ksea$ES)

# Extract time point from comparison
sig_ksea$timepoint <- str_match(string = sig_ksea$comparison, pattern = "PIC_(\\S+)-Control")[,2]


# Cast to wide format on logp
sig_ksea_wide <- dcast(data = sig_ksea, formula = kinase + kinase_gene ~ comparison, value.var = "logp")

# Get a matrix of logp values
logp_mat <- as.matrix(sig_ksea_wide[, 3:ncol(sig_ksea_wide)])

# Calculate distances
d <- dist(logp_mat)

# Cluster
clustered <- hclust(d)


# Plot as clustered heatmap
ggplot(data = sig_ksea, aes(x = sig_ksea$timepoint, y = sig_ksea$kinase_gene, color = sig_ksea$logp)) + geom_point(aes(size = sig_ksea$ES)) + scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "-LogP") + scale_y_discrete(limits = sig_ksea_wide[clustered$order, "kinase_gene"]) + scale_x_discrete(limits = c("05m", "30m", "02h", "24h")) + theme_bw() + scale_size(name = "ES") + xlab(label = "Time Point") + ylab(label = "Kinase")
ggsave(filename = "20230322-AK10-ph-ksea-dotplot.pdf", width = 4, height = 8)




# Make line plots of substrates for regulated kinases


# Create an empty table to store all sub results for plotting all at once
sub_results_all <- data.frame(
  Kinase = character(),
  Protein = character(),
  TimePoint = character(),
  log2FC = double(),
  Substrate = character()
)

# Loop through sig kinases
for (i in 1:length(sig_kinases)) {

  # Get current kinase
  current_kinase <- sig_kinases[i]

  # Get the gene name
  current_kinase_gene <- getGeneName[current_kinase]

  # Check that the kinase has a gene name
  if (!(is.na(current_kinase_gene))) {

    # Get  substrates in results
    subs <- protmapper[which(protmapper$CTRL_GENE_NAME == current_kinase_gene), "PH_ACC"]

    # Get results for substrates
    sub_results <- results[which(results$Protein %in% subs),]

    # Get any substrates with NA or Inf values to exclude
    exclude <- unique(sub_results[which(!(is.finite(sub_results$log2FC))), "Protein"])

    # Filter out excludes
    sub_results <- sub_results[which(!(sub_results$Protein %in% exclude)),]

    # Make a time point column that is easy to order for line plots (fix in illustrator later)
    sub_results[which(sub_results$Label == "PIC_05m-Control"), "TimePoint"] <- "1_05m"
    sub_results[which(sub_results$Label == "PIC_30m-Control"), "TimePoint"] <- "2_30m"
    sub_results[which(sub_results$Label == "PIC_02h-Control"), "TimePoint"] <- "3_02h"
    sub_results[which(sub_results$Label == "PIC_24h-Control"), "TimePoint"] <- "4_24h"

    # Add kinase gene name to sub rest
    sub_results$Kinase <- current_kinase_gene

    # Add Substrate column to indicate these are substrates
    sub_results$Substrate <- "Substrate"

    # Append sub_results data to sub_results_all
    sub_results_all <- rbind(sub_results[, c("Kinase", "Protein", "TimePoint", "log2FC", "Substrate")], sub_results_all)

    # Get average plots per time point
    sub_results_avg <- aggregate(log2FC ~ TimePoint, sub_results, mean)

    # Add a protein column
    sub_results_avg$Protein <- "Average"

    # Add Substrate column to indicate this is the average
    sub_results_avg$Substrate <- "Average"

    # Add kinase column
    sub_results_avg$Kinase <- current_kinase_gene

    # Append to sub_results_all
    sub_results_all <- rbind(sub_results_avg, sub_results_all)

    # Plot a line plot
    ggplot(data = sub_results, aes(x = sub_results$TimePoint, y = sub_results$log2FC)) + geom_line(group = sub_results$Protein) + theme_bw() + xlab(label = "Time Point") + ylab(label = "Log2FC") + ggtitle(label = paste(current_kinase_gene, "substrates")) + geom_line(data = sub_results_avg, group = sub_results_avg$Protein, aes(x = sub_results_avg$TimePoint, y = sub_results_avg$log2FC), color = "red", size = 2) + scale_x_discrete(expand = c(0,0))
    ggsave(filename = paste("./sigkinase_substrate_plots/", current_kinase_gene, ".pdf", sep = ""), width = 3, height = 3)
  }

}



# Make substrate plots all at once (this isn't working)
#ggplot(data = sub_results_all, aes(x = sub_results_all$TimePoint, y = sub_results_all$log2FC)) + geom_line(group = sub_results_all$Protein) + facet_wrap(Kinase ~ .)






# Plot MK02 185/187 ph over time

# Get results for MK01 185/187 ph
results_mk01 <- results[which(results$Protein == "P28482_ph185_ph187"),]

# Fix time points for order in line plotting
results_mk01[which(results_mk01$Label == "PIC_05m-Control"), "TimePoint"] <- "1_05m"
results_mk01[which(results_mk01$Label == "PIC_30m-Control"), "TimePoint"] <- "2_30m"
results_mk01[which(results_mk01$Label == "PIC_02h-Control"), "TimePoint"] <- "3_02h"
results_mk01[which(results_mk01$Label == "PIC_24h-Control"), "TimePoint"] <- "4_24h"

# Plot
ggplot(data = results_mk01, aes(x = results_mk01$TimePoint, y = results_mk01$log2FC)) + geom_line(group = results_mk01$Protein) + theme_bw() + scale_x_discrete(expand = c(0,0)) + xlab(label = "Time Point") + ylab(label = "Log2FC") + ggtitle(label = "MAPK1/ERK2 T185/Y187")
ggsave(filename = "MAPK1_T185_Y187_L2FC.pdf", width = 3, height = 3)



# Plot PLK1 210  ph over time

# Get results for MK01 185/187 ph
results_plk1 <- results[which(results$Protein == "P53350_ph210"),]

# Fix time points for order in line plotting
results_plk1[which(results_plk1$Label == "PIC_05m-Control"), "TimePoint"] <- "1_05m"
results_plk1[which(results_plk1$Label == "PIC_30m-Control"), "TimePoint"] <- "2_30m"
results_plk1[which(results_plk1$Label == "PIC_02h-Control"), "TimePoint"] <- "3_02h"
results_plk1[which(results_plk1$Label == "PIC_24h-Control"), "TimePoint"] <- "4_24h"

# Plot
ggplot(data = results_plk1, aes(x = results_plk1$TimePoint, y = results_plk1$log2FC)) + geom_line(group = results_plk1$Protein) + theme_bw() + scale_x_discrete(expand = c(0,0)) + xlab(label = "Time Point") + ylab(label = "Log2FC") + ggtitle(label = "PLK1 T210 Phos.")
ggsave(filename = "PLK1_T210_L2FC.pdf", width = 3, height = 3)
