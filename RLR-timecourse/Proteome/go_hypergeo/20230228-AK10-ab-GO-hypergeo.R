library(ggplot2)
library(reshape2)
library(stringr)
library(ggrepel)



# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital//Data/JohnsonLab/AndrewKurland/PolyIC-Timecourse/20230227-AK10-analysis/20230227-AK10-ab/go_hypergeo/")

# Set prefix for naming output files
prefix <- "20230321-AK10-ab-GO"

# Set DB file names
dbfiles <- c("~/Box Sync/DB/SwissProt.Human.2016.01.11.fasta")

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



# Read in results file
results <- read.delim(file = "../msstats/20230227-AK10-results-ann.txt", quote = "", stringsAsFactors = FALSE)


# Read in subquant file
subquant <- read.delim(file = "../msstats/20230227-AK10-subquant.txt", quote = "", stringsAsFactors = FALSE)


# Melt down subquant
subquant_melted <- melt(data = subquant, id.vars = "Protein", variable.name = "Sample", value.name = "Log2Intensity")

# Split sample into condition and replicate
subquant_melted$Condition <- str_match(string = subquant_melted$Sample, pattern = "(\\S+)_(\\d+)")[,2]
subquant_melted$BioReplicate <- str_match(string = subquant_melted$Sample, pattern = "(\\S+)_(\\d+)")[,3]

# Replace NA values with 0s
subquant_melted[which(is.na(subquant_melted$Log2Intensity)), "Log2Intensity"] <- 0

# Aggregate by getting max/min of Log2Intensity for each protein in each condition
subquant_agg_max <- aggregate(Log2Intensity ~ Protein + Condition, data = subquant_melted, FUN = max)
subquant_agg_min <- aggregate(Log2Intensity ~ Protein + Condition, data = subquant_melted, FUN = min)

# Make lookup vectors to get max/min Log2Intensity of each protein in each condition
getSubquantMax <- subquant_agg_max$Log2Intensity
names(getSubquantMax) <- paste(subquant_agg_max$Protein, subquant_agg_max$Condition)
getSubquantMin <- subquant_agg_min$Log2Intensity
names(getSubquantMin) <- paste(subquant_agg_min$Protein, subquant_agg_min$Condition)



# Set cutoff values for signifiance
l2fc_max <- 1
l2fc_min <- -1
pval_max <- 0.05
adjpval_max <- 0.05



# How to treat Inf values
# "ignore" removes them, "include" includes them, "strict" only includes them if they are present in all conditions of the one condition being compared and none of the other
infinclude <- "ignore"



# Read in GO term definitions
go_def <-
  read.delim(file = "~/Box Sync/DB/GO/20210218-GO-term-definitions.txt",
             quote = "",
             header = FALSE)


# Created named vector to look up GO definitions
getGOdefinition <- go_def$V3
names(getGOdefinition) <- go_def$V1

# Read in GO annotations by UniProt ID
go_ann <-
  read.delim(file = "~/Box Sync/DB/GO/20210218-UniProt-GO-annotations.txt", header = FALSE)



# Create data frame to store enrichment results (pvalues)
go_enrichment_table <- data.frame(
  comparison = character(),
  direction = character(),
  go_term = character(),
  go_definition = character(),
  x = integer(),
  m = integer(),
  n = integer(),
  k = integer(),
  pvalue = double(),
  logp = double(),
  proteins = character(),
  protein_names = character()
)




# Get list of comparisons in results file
comparisons <- unique(results$Label)


# Loop through comparisons one at a time and test GO enrichments
for (i in 1:length(comparisons)) {

  # Get current comparison
  currentcomparison <- comparisons[i]

  # Split comparison by hyphens into conditions
  current_condition1 <- unlist(str_split(currentcomparison, "-"))[1]
  current_condition2 <- unlist(str_split(currentcomparison, "-"))[2]

  print(currentcomparison)


  # Get list of upregulated proteins for this comparison
  upreg_proteins <-
    unique(results[which((results$Label == currentcomparison) &
                           (results$log2FC > l2fc_max) &
                           (results$pval < pval_max) &
                           (results$adj.pvalue < adjpval_max) &
                           (is.finite(results$log2FC)) &
                           (is.finite(results$pvalue))
    ), "Protein.Accession"])


  # Get proteins with positive Inf l2fc values
  posinf_proteins <- unique(results[which(
    (results$Label == currentcomparison) &
      (results$log2FC == Inf)), "Protein.Accession"])



  # Get proteins with Inf l2fc values present in all condition1 and none of condition2
  posinf_proteins_strict <- subquant[which(
    (subquant$Protein %in% posinf_proteins) &
      (getSubquantMax[paste(subquant$Protein, current_condition2)] == 0) &
      (getSubquantMin[paste(subquant$Protein, current_condition1)] > 0)), "Protein.Accession"]


  # Add +Inf proteins to upreg proteins depending on how infinclude is set above
  if (infinclude == "include") {
    upreg_proteins <- append(upreg_proteins, posinf_proteins)
  } else if (infinclude == "strict") {
    upreg_proteins <- append(upreg_proteins, posinf_proteins_strict)
  }


  # Get list of downregulated proteins for this comparison
  downreg_proteins <-
    unique(results[which((results$Label == currentcomparison) &
                           (results$log2FC < l2fc_min) &
                           (results$pvalue < pval_max) &
                           (results$adj.pvalue < adjpval_max) &
                           (is.finite(results$log2FC)) &
                           (is.finite(results$pvalue))
    ), "Protein.Accession"])



  # Get proteins with negative Inf l2fc values
  neginf_proteins <- unique(results[which(
    (results$Label == currentcomparison) &
      (results$log2FC == -Inf)), "Protein.Accession"])

  # Get proteins with -Inf l2fc values present in all condition2 and none of condition1
  neginf_proteins_strict <- subquant[which(
    (subquant$Protein %in% neginf_proteins) &
      (getSubquantMax[paste(subquant$Protein, current_condition1)] == 0) &
      (getSubquantMin[paste(subquant$Protein, current_condition2)] > 0)), "Protein.Accession"]

  # Add -Inf proteins to dwonreg proteins depending on how infinclude is set above
  if (infinclude == "include") {
    downreg_proteins <- unique(append(downreg_proteins, neginf_proteins))
  } else if (infinclude == "strict") {
    downreg_proteins <- unique(append(downreg_proteins, neginf_proteins_strict))
  }

  # Get list of all proteins detected for this comparison
  all_proteins <- unique(append(results[which(is.finite(results$log2FC)), "Protein.Accession"], append(upreg_proteins, downreg_proteins)))


  # Get GO terms containing up or down reuglated proteins to test
  go_terms <-
    unique(go_ann[which((go_ann$V2 %in% upreg_proteins) |
                          (go_ann$V2 %in% downreg_proteins)), 5])

  # Loop through each GO term and test for enrichment in up/down regulated proteins
  if (length(go_terms) > 0) {
    for (j in 1:length(go_terms)) {

      # Get the current GO term
      currentterm <- go_terms[j]

      # Get the definition for the current term
      current_def <- unname(getGOdefinition[as.character(currentterm)])


      # Get list of proteins in the current term
      go_proteins <-
        go_ann[which(go_ann$V5 == currentterm), "V2"]


      # Only continue if # proteins in this go term is < 100
      if (length(go_proteins) < 100) {

        # Perform hypergeomteric test on upreg proteins first

        # Get number of white balls drawn from the urn (upreg in GO term)
        x <- length(intersect(upreg_proteins, go_proteins))

        # Only continue if x >= 2
        if (x >= 2) {

          # Get number of white balls in the urn (number of protein with current term in full dataset)
          m <- length(intersect(all_proteins, go_proteins))

          # Get the number of black balls in the urn
          n <- length(all_proteins) - m

          # Get the number of balls drawn from the urn (number of upreg proteins)
          k <- length(upreg_proteins)

          # Calculate p-value by hypergeometric test
          pvalue <- dhyper(x = x,
                           m = m,
                           n = n,
                           k = k)


          # Calculate log 10 pvalue
          logp <- -1 * log(pvalue) / log(10)

          # Get list of protein accession upregulated for this GO term
          upreg_currentterm <-
            paste(
              intersect(upreg_proteins, go_proteins),
              sep = "",
              collapse = ";"
            )


          # GEt protein names for these proteins
          protein_names_currentterm <-
            paste(
              unname(getProteinName[intersect(upreg_proteins, go_proteins)]),
              sep = "",
              collapse = ";"
            )

          # Bind results of enrichment test to the pvalue table
          go_enrichment_table <-
            rbind(
              data.frame(
                comparison = currentcomparison,
                direction = "UP",
                go_term = currentterm,
                go_definition = current_def,
                x = x,
                m = m,
                n = n,
                k = k,
                pvalue = pvalue,
                logp = logp,
                proteins = upreg_currentterm,
                protein_names = protein_names_currentterm
              ),
              go_enrichment_table
            )

        }


        # Now do the same for downreuglated proteins

        # Perform hypergeomteric test on downreg proteins
        # Get number of which balls drawn from the urn (upreg in GO term)
        x <- length(intersect(downreg_proteins, go_proteins))



        # Only continue if x >= 2
        if (x >= 2) {

          # Get number of white balls in the urb (number of protein with current term in full dataset)
          m <- length(intersect(all_proteins, go_proteins))

          # Get the number of black balls in the urn
          n <- length(all_proteins) - m

          # Get the number of balls drawn from the urn (number of upreg proteins)
          k <- length(downreg_proteins)

          # Calculate p-value by hypergeometric test
          pvalue <- dhyper(x = x,
                           m = m,
                           n = n,
                           k = k)

          # Calculate log 10 pvalue
          logp <- -1 * log(pvalue) / log(10)

          # Get list of protein accession upregulated for this GO term
          downreg_currentterm <-
            paste(
              intersect(downreg_proteins, go_proteins),
              sep = "",
              collapse = ";"
            )


          # GEt protein names for these proteins
          protein_names_currentterm <-
            paste(
              unname(getProteinName[intersect(downreg_proteins, go_proteins)]),
              sep = "",
              collapse = ";"
            )

          # Bind results of enrichment test to the pvalue table
          go_enrichment_table <-
            rbind(
              data.frame(
                comparison = currentcomparison,
                direction = "DOWN",
                go_term = currentterm,
                go_definition = current_def,
                x = x,
                m = m,
                n = n,
                k = k,
                pvalue = pvalue,
                logp = logp,
                proteins = downreg_currentterm,
                protein_names = protein_names_currentterm
              ),
              go_enrichment_table
            )
        }
      }
    }
  }
}



# Adjust pvalues for multiple testing
go_enrichment_table$adj.pvalue <-
  p.adjust(p = go_enrichment_table$pvalue, method = "fdr")

# Write GO enrichment table to file
write.table(
  x = go_enrichment_table,
  file = "20230228-AK10-ab-GO-hypergeo-adjpval05-infIgnore.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# Filter for significant terms by adj.pval < 0.05
sig_terms <- unique(go_enrichment_table[which(go_enrichment_table$adj.pvalue < 0.05), "go_term"])

# Filter go table for sig terms
go_sig <- go_enrichment_table[which(go_enrichment_table$go_term %in% sig_terms),]


# Create an empty list to store top terms
top_terms <- list()


# Loop through comparisons and get top 10 terms from each time point
for (i in 1:length(comparisons)) {

  # Get the current comparison
  current_comparison <- comparisons[i]

  # Filter go_sig for current comparison
  current_go <- go_sig[which(go_sig$comparison == current_comparison),]

  # Order the table by pvalue
  current_go <- current_go[order(current_go$pvalue),]

  # If there are fewer than 10 lines then take all
  if (nrow(current_go) < 10) {

    # Get the current terms to add to the top_terms list
    current_terms <- current_go$go_term

    # Append these to top_terms list
    top_terms <- append(top_terms, current_terms)
  } else {

    # Otherwise add only the top 10 terms
    current_terms <- current_go[1:10, "go_term"]

    # Append these to top_terms list
    top_terms <- append(top_terms, current_terms)
  }

}

# Get unique top terms
top_terms <- unique(unlist(top_terms))


# Filter go table for top terms
top_go <- go_enrichment_table[which(go_enrichment_table$go_term %in% top_terms),]

# Remove IFN gamma term downreg from 24h
top_go <- top_go[which(!((top_go$go_term == "GO:0060333") & (top_go$direction == "DOWN"))),]

# Flip the sign on logp for downreg term
top_go[which(top_go$direction == "DOWN"), "logp"] <- -1 * top_go[which(top_go$direction == "DOWN"), "logp"]


# Cast to wide format on logp for clustering
top_go_wide <- dcast(data = top_go, formula = go_definition + go_term ~ comparison, value.var = "logp")

# Replace NAs with zeros
top_go_wide[is.na(top_go_wide)] <- 0

# Convert to a matrix
logp_mat <- as.matrix(top_go_wide[,3:ncol(top_go_wide)])

# Calculate distances
d <- dist(logp_mat)

# Cluster
clustered <- hclust(d)

# Get the order of go_definitions for plotting
go_def_order <- top_go_wide[clustered$order, "go_definition"]


# Split into up and downregulated terms
#top_go_up <- top_go[which(top_go$direction == "UP"),]
#top_go_down <- top_go[which(top_go$direction == "DOWN"),]


# Cast top_go_up to wide format on LogP for clustering
#top_go_up_wide <- dcast(data = top_go_up, formula = go_definition + go_term ~ comparison, value.var = "logp")

# Replace NAs with zeros
#top_go_up_wide[is.na(top_go_up_wide)] <- 0

# Get matrix of logps for clustering
#logp_mat <- as.matrix(top_go_up_wide[,3:ncol(top_go_up_wide)])

# Calculate distances
#d <- dist(logp_mat)

# Cluster
#clustered <- hclust(d)

# Get the order of go_definitions for plotting
#go_def_order <- top_go_up_wide[clustered$order, "go_definition"]


# PLot upreg terms as dot plot
ggplot(data = top_go, aes(x = top_go$comparison, y = top_go$go_definition)) + geom_point(aes(size = top_go$x, color = top_go$logp)) + scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "-LogP") + theme_bw() + scale_y_discrete(limits = go_def_order) + scale_x_discrete(limits = c("PIC_02h-Control", "PIC_24h-Control")) + theme(axis.text.x = element_text(angle = 90)) + xlab(label = NULL) + ylab(label = NULL) + scale_size(name = "Overlap") + ggtitle(label = "Top GO terms enriched with |L2FC| > 1 & adj.pval < 0.05")
ggsave(filename = "20230228-AK10-ab-GO-hypergeo-dotplot.pdf", width = 6, height = 6)




















# Cast top_go_up to wide format on LogP for clustering
top_go_down_wide <- dcast(data = top_go_down, formula = go_definition + go_term ~ comparison, value.var = "logp")

# Replace NAs with zeros
top_go_down_wide[is.na(top_go_down_wide)] <- 0

# Get matrix of logps for clustering
logp_mat <- as.matrix(top_go_down_wide[,3:ncol(top_go_down_wide)])

# Calculate distances
d <- dist(logp_mat)

# Cluster
clustered <- hclust(d)

# Get the order of go_definitions for plotting
go_def_order <- top_go_down_wide[clustered$order, "go_definition"]


# PLot downreg terms as dot plot
ggplot(data = top_go_down, aes(x = top_go_down$comparison, y = top_go_down$go_definition)) + geom_point(aes(size = top_go_down$x, color = top_go_down$logp)) + scale_color_gradient(low = "white", high = "blue", name = "-LogP") + theme_bw() + scale_y_discrete(limits = go_def_order) + scale_x_discrete(limits = c("PIC_05m-Control", "PIC_30m-Control", "PIC_02h-Control", "PIC_24h-Control")) + theme(axis.text.x = element_text(angle = 90)) + xlab(label = NULL) + ylab(label = NULL) + scale_size(name = "Overlap") + ggtitle(label = "GO terms enriched with L2FC < -1 & adj.pval < 0.05")
ggsave(filename = "20221116-AK10-ab-GO-hypergeo-downreg-dotplot.pdf", width = 8, height = 6)




# Plot protein abundance plots for each top term
for (i in 1:length(top_go_up)) {

  # Get current term and definiton
  current_term <- top_go_up[i, "go_term"]
  current_def <- top_go_up[i, "go_definition"]

  # Get proteins in current terms
  current_proteins <- top_go_up[i, "proteins"]

  # Split by semicolons
  current_proteins_list <- unlist(str_split(string = current_proteins, pattern = ";"))


  # Extract current proteins from results
  current_results <- results[which(results$Protein %in% current_proteins_list),]

  # Replace Inf/-Inf values with NA
  current_results[which(current_results$log2FC == "-Inf"), "log2FC"] <- NA
  current_results[which(current_results$log2FC == "Inf"), "log2FC"] <- NA

  # Plot as line plots
  ggplot(data = current_results, aes(x = current_results$Label, y = current_results$log2FC)) + geom_line(aes(color = current_results$Protein.Name, group = current_results$Protein.Name)) + scale_x_discrete(limits = c("PIC_05m-Control", "PIC_30m-Control", "PIC_02h-Control", "PIC_03h-Control", "PIC_04h-Control", "PIC_08h-Control", "PIC_24h-Control")) + theme_bw() + scale_color_discrete(name = NULL) + geom_point(aes(color = current_results$Protein.Name)) + ggtitle(label = current_def) + xlab(label = "Time Point") + ylab(label = "Log2FC") + theme(axis.text.x = element_text(angle = 90))
  ggsave(filename = paste("./Plots/", current_def, ".pdf", sep = ""), width = 6, height = 6)


}



# Plot protein abundance plots for each top term
for (i in 1:length(top_go_down)) {

  # Get current term and definiton
  current_term <- top_go_down[i, "go_term"]
  current_def <- top_go_down[i, "go_definition"]

  # Get proteins in current terms
  current_proteins <- top_go_down[i, "proteins"]

  # Split by semicolons
  current_proteins_list <- unlist(str_split(string = current_proteins, pattern = ";"))


  # Extract current proteins from results
  current_results <- results[which(results$Protein %in% current_proteins_list),]

  # Replace Inf/-Inf values with NA
  current_results[which(current_results$log2FC == "-Inf"), "log2FC"] <- NA
  current_results[which(current_results$log2FC == "Inf"), "log2FC"] <- NA

  # Plot as line plots
  ggplot(data = current_results, aes(x = current_results$Label, y = current_results$log2FC)) + geom_line(aes(color = current_results$Protein.Name, group = current_results$Protein.Name)) + scale_x_discrete(limits = c("PIC_05m-Control", "PIC_30m-Control", "PIC_02h-Control", "PIC_03h-Control", "PIC_04h-Control", "PIC_08h-Control", "PIC_24h-Control")) + theme_bw() + scale_color_discrete(name = NULL) + geom_point(aes(color = current_results$Protein.Name)) + ggtitle(label = current_def) + xlab(label = "Time Point") + ylab(label = "Log2FC") + theme(axis.text.x = element_text(angle = 90))
  ggsave(filename = paste("./Plots/", current_def, ".pdf", sep = ""), width = 6, height = 6)


}





# Get sig GO terms at 24h for bar plotting
sig_go_terms <- go_enrichment_table[which(
  (go_enrichment_table$adj.pvalue < 0.05) &
    (go_enrichment_table$comparison == "PIC_24h-Control")), "go_term"]

# Filter for sig go terms
go_sig <- go_enrichment_table[which(go_enrichment_table$go_term %in% sig_go_terms),]


# Get upreg terms
go_sig_up <- go_sig[which(go_sig$direction == "UP"),]

# order by pvalue
go_sig_up <- go_sig_up[order(go_sig_up$pvalue),]

# Get the top 10 terms
top10_up_terms <- go_sig_up[1:10, "go_term"]

# Get up terms from go enrichment table (24 h only)
top10_up_table <- go_enrichment_table[which(
  (go_enrichment_table$go_term %in% top10_up_terms) &
    (go_enrichment_table$direction == "UP") &
    (go_enrichment_table$comparison == "PIC_24h-Control")),]


# Get downreg terms
go_sig_down <- go_sig[which(go_sig$direction == "DOWN"),]

# order by pvalue
go_sig_down <- go_sig_down[order(go_sig_down$pvalue),]

# Get the top 10 terms
top10_down_terms <- go_sig_down[1:10, "go_term"]


# Get up terms from go enrichment table (24 h only)
top10_down_table <- go_enrichment_table[which(
  (go_enrichment_table$go_term %in% top10_down_terms) &
    (go_enrichment_table$direction == "DOWN") &
    (go_enrichment_table$comparison == "PIC_24h-Control")),]

# Combine top10 table
top10_table <- rbind(top10_up_table, top10_down_table)

# Flip the sign on LogP for downreg terms
top10_table[which(top10_table$direction == "DOWN"), "logp"] <- -1 * top10_table[which(top10_table$direction == "DOWN"), "logp"]

# Set colors for bar plots
colors <- c("UP" = "red", "DOWN" = "blue")

# Plot as a horizontal bar plot
ggplot(data = top10_table, aes(x = top10_table$go_definition, y = top10_table$logp)) + geom_bar(stat = "identity", aes(fill = top10_table$direction)) + coord_flip() + scale_x_discrete(limits = top10_table[order(top10_table$logp), "go_definition"]) + theme_bw() + scale_fill_manual(values = colors, name = "Direction") + ylab(label = "-Log10(p-value)") + xlab(label = "")
ggsave(filename = paste(prefix, "-GO-24h-barplot.pdf", sep = ""), width = 6, height = 4)
