library(ggplot2)
library(reshape2)
library(stringr)
library(gplots)

# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital/Data/JohnsonLab/AndrewKurland/PolyIC-siRNA/go-hypergeo/")


# Set prefix
prefix <- "20230329-AKsiRNA-GO-hypergeo"

# Set DB file names
dbfiles <- c("~/OneDrive - The Mount Sinai Hospital/DB/20191010.SwissProt.Hsapiens.fasta")

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
results <- read.delim(file = "../msstats/20230327-AKsiRNA-results-ann.txt", quote = "", stringsAsFactors = FALSE)


# Read in subquant file
subquant <- read.delim(file = "../msstats/20230327-AKsiRNA-subquant.txt", quote = "", stringsAsFactors = FALSE)



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
  file = paste(prefix, "-adjpval05-infIgnore.txt", sep = ""),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)




# Make a list of comparisons of interest
comparisons_of_interest <- c("siATR-siControl", "siATM-siControl", "siAMPKA1-siControl")

# Make a column combining term and direction
go_enrichment_table$term_direction <- paste(go_enrichment_table$go_term, go_enrichment_table$direction)

# Make a list for top terms
top_terms <- list()


# Loop through comparisons of interest
for (i in 1:length(comparisons_of_interest)) {

  # Get current comparison
  current_comparison <- comparisons_of_interest[i]

  # Get go results for current comaprison
  current_go <- go_enrichment_table[which(go_enrichment_table$comparison == current_comparison),]

  # Filter by adjpval < 0.05
  current_go <- current_go[which(current_go$adj.pvalue < 0.05),]

  # Sort by pvalue
  current_go <- current_go[order(current_go$pvalue, decreasing = FALSE),]

  # Check if there are at least 10 terms
  if (nrow(current_go) > 10) {

    # Get top 10 terms
    top10_terms <- current_go[1:10, "term_direction"]

    # Append to top_terms list
    top_terms <- append(top_terms, top10_terms)

  } else {

    # Otherwise get all the terms
    top10_terms <- current_go[, "term_direction"]

    # Append to top terms list
    top_terms <- append(top_terms, top10_terms)

  }


}



# Filter for unique top terms
top_terms <- unique(top_terms)



# Filter go results for comparisons of interest and top terms
go_filtered <- go_enrichment_table[which(
  (go_enrichment_table$comparison %in% comparisons_of_interest) &
    (go_enrichment_table$term_direction %in% top_terms)), ]


# Flip the sign on logp for downreg terms
go_filtered[which(go_filtered$direction == "DOWN"), "logp"] <- -1 * go_filtered[which(go_filtered$direction == "DOWN"), "logp"]

# Cast to wide format on logp
go_filtered_wide <- dcast(data = go_filtered, formula = go_term + go_definition ~ comparison, value.var = "logp")

# Get numeric data as matrix
logp_mat <- as.matrix(go_filtered_wide[, 3:ncol(go_filtered_wide)])

# Replace NAs with 0s
logp_mat[is.na(logp_mat)] <- 0

# Calculate distances
d <- dist(logp_mat)

# Cluster
clustered <- hclust(d)


# Plot as clustered dot plot
ggplot(data = go_filtered, aes(x = go_filtered$comparison, y = go_filtered$go_definition)) + geom_point(aes(color = go_filtered$logp, size = go_filtered$x)) + scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "-LogP") + scale_y_discrete(limits = go_filtered_wide[clustered$order, "go_definition"]) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + xlab(label = "Comparison") + ylab(label = "GO definition") + scale_size(name = "Overlap")
ggsave(filename = paste(prefix, "-dotplot.pdf", sep = ""), width = 6, height = 7)
