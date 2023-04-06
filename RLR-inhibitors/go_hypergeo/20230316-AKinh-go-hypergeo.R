library(ggplot2)
library(reshape2)
library(stringr)
library(gplots)



# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital/Data/JohnsonLab/AndrewKurland/PolyIC-Inhibitor-Screen/20230316-AKinh-go-hypergeo/")

# Set prefix
prefix <- "20230316-AKinh-GO-hypergeo"

# Set DB file names
dbfiles <- c("~/OneDrive - The Mount Sinai Hospital/DB/SwissProt.Human.2016.01.11.fasta")

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
results <- read.delim(file = "../20230315-AKinh-msstats/20230315-AK-inh-ab-results-ann.txt", quote = "", stringsAsFactors = FALSE)


# Read in subquant file
subquant <- read.delim(file = "../20230315-AKinh-msstats/20230315-AK-inh-ab-subquant.txt", quote = "", stringsAsFactors = FALSE)


# Melt down subquant
subquant_melted <- melt(data = subquant, id.vars = "Protein", variable.name = "Sample", value.name = "Log2Intensity")

# Split sample into condition and replicate
subquant_melted$Condition <- str_match(string = subquant_melted$Sample, pattern = "(\\S+)\\s*_\\d+")[,2]
subquant_melted$BioReplicate <- str_match(string = subquant_melted$Sample, pattern = "(\\S+)\\s*_(\\d+)")[,3]


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
  read.delim(file = "~/OneDrive - The Mount Sinai Hospital//DB/GO/20210218-GO-term-definitions.txt",
             quote = "",
             header = FALSE)


# Created named vector to look up GO definitions
getGOdefinition <- go_def$V3
names(getGOdefinition) <- go_def$V1

# Read in GO annotations by UniProt ID
go_ann <-
  read.delim(file = "~/OneDrive - The Mount Sinai Hospital//DB/GO/20210218-UniProt-GO-annotations.txt", header = FALSE)


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
  file = "20230316-AKinh-GO-hypergeo-adjpval05-infIgnore.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)



# Filter for significant terms by adj.pval < 0.01 and x > 4
sig_terms <- unique(go_enrichment_table[which(go_enrichment_table$adj.pvalue < 0.01), "go_term"])

# Filter go table for sig terms
go_sig <- go_enrichment_table[which(go_enrichment_table$go_term %in% sig_terms),]

# Drop the DMSO_Control comparison for this, it's alreaddy in Figure 1
go_sig <- go_sig[which(!(go_sig$comparison == "DMSO_Control")),]

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

# Drop DMSO_Control comparison again
top_go <- top_go[which(!(top_go$comparison == "DMSO_Control")),]

# Cast to wide format on logp
top_go_wide <- dcast(top_go, formula = go_definition + go_term + direction ~ comparison, value.var = "logp")

# Get downreg rows
downreg_rows <- which(top_go_wide$direction == "DOWN")

# Loop through downreg row
for (i in 1:length(downreg_rows)) {

  # Get current row
  current_row <- downreg_rows[i]

  # Loop through columns of current row
  for (j in 4:ncol(top_go_wide)) {
    top_go_wide[current_row, j] <- -1 * top_go_wide[current_row, j]
  }

}

# Get numeric data for heatmap plotting
logp_mat <- as.matrix(top_go_wide[, 4:ncol(top_go_wide)])
row.names(logp_mat) <- top_go_wide$go_definition

# Replace NAs with zeros
logp_mat[is.na(logp_mat)] <- 0

# Set color scale from red to white to blue from 0 to 1 in 0.1 increments
breaks <- seq(-1, 1, 0.1)
colorpalette <- colorRampPalette(c("blue", "white", "red"))

# Plot heatmap of correlations and save to pdf
pdf(paste(prefix, "-GO-hypergeo-heatmap-large.pdf", collapse = '', sep = ''), width = 10, height = 20)
heatmap.2(x=logp_mat, col = colorpalette(100), trace="none", density.info = "none", key.xlab = "-LogP", margins = c(10,20))
dev.off()
