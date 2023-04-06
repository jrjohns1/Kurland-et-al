library(reshape2)
library(stringr)
library(seqinr)
library(ggseqlogo)

# Set working directory
setwd("~/Box Sync/Data/JohnsonLab/AndrewKurland/PolyIC-Timecourse/20230227-AK10-analysis/20230227-AK10-ph/motif-x/")


# Read in results
results <- read.delim(file = "../msstats/20230227-AK10-ph-results-ann.txt")

# Set db files
dbfiles <- "~/Box Sync/DB/SwissProt.Human.2016.01.11.fasta"

numdb <- length(dbfiles)

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

# Set a flag to zero for when to read in sequences
grabsequence <- 0

# Declare lists to store sequences and accessions
sequencelist <- list()
accessionlist <- list()

# Loop through databases line by line
for (i in 1:length(dblines)) {
  # If the line contains a ">" then this is a description line
  if(grepl(pattern = ">", x = dblines[i])) {
    # If this is not the first line then append the sequence and accession to their respective lists
    if (i > 1) {
      sequencelist <- append(sequencelist, sequence)
      accessionlist <- append(accessionlist, accession)
    }
    # Extract accession based on the pattern: ">sp|P12345|ASDF_HUMAN blah blah", where P12345 is accession
    accession <- str_match(string = dblines[i], pattern = ">\\S\\S\\|([^\\|]*)\\|.*")[2]

    # Set grab sequence flag to 1
    grabsequence <- 1
  } else {
    if (grabsequence == 1) {
      # If this is not a header line and grabsequence flag is set then read line into sequence
      sequence <- dblines[i]
      # And set the grabsequence flag back to zero
      grabsequence <- 0
    } else {
      # Otherwise append the sequence to the existing sequqence
      sequence <- paste(sequence, dblines[i], sep = '', collapse = '')
    }
  }
}

# Grab the last ones in the file
sequencelist <- append(sequencelist, sequence)
accessionlist <- append(accessionlist, accession)

# Store sequences in a named vector for quick retrieval by accession
getSequence <- sequencelist
names(getSequence) <- accessionlist












# Get phaccs
phaccs <- as.data.frame(unique(results$Protein))
colnames(phaccs) <- "Accession"

###################

# Get flanking sequences


# Loop through list of phaccs
for (i in 1:nrow(phaccs)) {

  # Check if this is an unambiguous ID
  if (str_detect(string = phaccs$Accession[i], ";")) {

    # Split by semi colons
    tmp <- str_split(string = phaccs$Accession[i], pattern = ";")[[1]][1]

    # Grab the first one to use to find flankseq
    phacc_split <- str_split(string = tmp, pattern = "_ph")
  } else {
    # Otherwise just split ph accession on _ph characters
    phacc_split <- str_split(string = phaccs$Accession[i], pattern = "_ph")

  }


  # Get protein accession
  protein.accession <- phacc_split[[1]][1]

  # Loop through sites
  for (j in 2:length(phacc_split[[1]])) {

    # Get sequence for this protein
    sequence <- getSequence[protein.accession]

    # Get the length of the sequence
    seqlength <- str_length(sequence)

    # Get phosphosite position
    position <- as.numeric(phacc_split[[1]][j])

    # Calculate start and end of positions in string to grab
    start <- position - 7
    end <- position + 7

    # Check that there are at least 7 AAs before and after position to grab
    if (start < 1) {
      start <- 1
    }
    if (end > seqlength) {
      end = seqlength
    }

    # First grab the sequence up to and including the phosphosite
    flankseq <- str_sub(string = sequence, start = start, end = position)

    # Add a * to denote the phospho position
    flankseq <- paste(flankseq, "*", sep = "", collapse = "")

    # Extract sequence
    flankseq <- paste(flankseq, str_sub(string = sequence, start = position + 1, end = end), sep = "", collapse = "")



    # If this is a multiply phosphorylated site then concatenate by | characters
    if (j > 2) {
      flankseq_out <- paste(flankseq_out, flankseq, sep = "|")
    } else {
      flankseq_out <- flankseq
    }



  }
  # Add flankseq to phacc dataframe
  phaccs[i, "Flankseq"] <- flankseq_out




}


# Create a named vector of flankseqs
getFlankseq <- phaccs$Flankseq
names(getFlankseq) <- phaccs$Accession



# Add Flankseq to results
results$Flankseq <- getFlankseq[results$Protein]

# Write to a file
write.table(x = results, file = "2020228-AK10-ph-results-flankseq.txt", sep = "\t", row.names = FALSE, quote = FALSE)




# Filter out flankseq with length other than 16
results <- results[which(str_length(results$Flankseq) == 16),]

# Remove *s in flankseq
results$Flankseq <- str_replace(string = results$Flankseq, pattern = "\\*", replacement = "")


# Set values for diff abundance
l2fc_min <- -1
l2fc_max <- 1
pval_max <- 0.05
adjpval_max <- 0.05

# Get a list of comparisons
comparisons <- unique(results$Label)


# Loop through comparisons and write out files for motif-x analysis
for (i in 1:length(comparisons)) {

  # Get the current comparison
  current_comparison <- comparisons[i]

  # Get list of upreg flankseq
  upreg_flankseq <- results[which(
    (results$Label == current_comparison) &
      (results$log2FC > l2fc_max) &
      (results$pvalue < pval_max) &
      (results$adj.pvalue < adjpval_max)), "Flankseq"]

  # Get names
  upreg_names <- results[which(
    (results$Label == current_comparison) &
      (results$log2FC > l2fc_max) &
      (results$pvalue < pval_max) &
      (results$adj.pvalue < adjpval_max)), "Protein"]

  # Write out as fasta file
  write.fasta(sequences = as.list(upreg_flankseq), names = upreg_names, file.out = paste(current_comparison, "-upreg.fasta", sep = ""))

  # Get list of downreg flankseq
  downreg_flankseq <- results[which(
    (results$Label == current_comparison) &
      (results$log2FC < l2fc_min) &
      (results$pvalue < pval_max) &
      (results$adj.pvalue < adjpval_max)), "Flankseq"]

  # Get names
  downreg_names <- results[which(
    (results$Label == current_comparison) &
      (results$log2FC < l2fc_min) &
      (results$pvalue < pval_max) &
      (results$adj.pvalue < adjpval_max)), "Protein"]

  # Write out as fasta file
  write.fasta(sequences = as.list(downreg_flankseq), names = downreg_names, file.out = paste(current_comparison, "-downreg.fasta", sep = ""))





  # Get list of background flankseq
  all_flankseq <- results[which(results$Label == current_comparison), "Flankseq"]

  # GEt all names
  all_names <- results[which(results$Label == current_comparison), "Protein"]

  # Write out as fasta file
  write.fasta(sequences = as.list(all_flankseq), names = all_names, file.out = paste(current_comparison, "-all.fasta", sep = ""))

}









# Read in 30m upreg motif-x output files
results_30m_upreg <- read.delim(file = "30m-upreg/momo.tsv", quote = "", stringsAsFactors = FALSE)

# Filter out comment lines
results_30m_upreg <- results_30m_upreg[which(str_length(results_30m_upreg$motif) > 0),]

# Add column with time point
results_30m_upreg$timepoint <- "30m"

# Add column with direction
results_30m_upreg$direction <- "upreg"



# Read in 30m downreg motif-x output files
results_30m_downreg <- read.delim(file = "30m-downreg/momo.tsv", quote = "", stringsAsFactors = FALSE)

# Filter out comment lines
results_30m_downreg <- results_30m_downreg[which(str_length(results_30m_downreg$motif) > 0),]

# Add column with time point
results_30m_downreg$timepoint <- "30m"

# Add column with direction
results_30m_downreg$direction <- "downreg"



# Read in 2h upreg motif-x output files
results_2h_upreg <- read.delim(file = "2h-upreg/momo.tsv", quote = "", stringsAsFactors = FALSE)

# Filter out comment lines
results_2h_upreg <- results_2h_upreg[which(str_length(results_2h_upreg$motif) > 0),]

# Add column with time point
results_2h_upreg$timepoint <- "2h"

# Add column with direction
results_2h_upreg$direction <- "upreg"


# Read in 2h downreg motif-x output files
results_2h_downreg <- read.delim(file = "2h-downreg/momo.tsv", quote = "", stringsAsFactors = FALSE)

# Filter out comment lines
results_2h_downreg <- results_2h_downreg[which(str_length(results_2h_downreg$motif) > 0),]

# Add column with time point
results_2h_downreg$timepoint <- "2h"

# Add column with direction
results_2h_downreg$direction <- "downreg"




# Read in 24h upreg motif-x output files
results_24h_upreg <- read.delim(file = "24h-upreg/momo.tsv", quote = "", stringsAsFactors = FALSE)

# Filter out comment lines
results_24h_upreg <- results_24h_upreg[which(str_length(results_24h_upreg$motif) > 0),]

# Add column with time point
results_24h_upreg$timepoint <- "24h"

# Add column with direction
results_24h_upreg$direction <- "upreg"


# Read in 24h downreg motif-x output files
results_24h_downreg <- read.delim(file = "24h-downreg/momo.tsv", quote = "", stringsAsFactors = FALSE)

# Filter out comment lines
results_24h_downreg <- results_24h_downreg[which(str_length(results_24h_downreg$motif) > 0),]

# Add column with time point
results_24h_downreg$timepoint <- "24h"

# Add column with direction
results_24h_downreg$direction <- "downreg"




# Combine all motif-x results together
results_motifx <- rbind(results_30m_upreg, results_30m_downreg, results_2h_upreg, results_2h_downreg, results_24h_upreg, results_24h_downreg)

# Save to a file
write.table(results, file  = "20230301-motif-x-results.txt", sep = "\t", quote = FALSE, row.names = FALSE)







# Make a list of motifs of interest
motifs_of_interest <- c("..........E....", ".........E.....", "........D......", "........D.E....", "........F......", "........L......", "........P......", "....R..........")


# Make a data frame to store enrichment results
motif_enrichment_table <- data.frame(
  comparison = character(),
  motif = character(),
  x = integer(),
  m = integer(),
  n = integer(),
  k = integer(),
  pvalue = numeric(),
  direction = character()

)



# Loop through comparisons
for (i in 1:length(comparisons)) {

  # Get current comparison
  current_comparison <- comparisons[i]

  # Loop through motifs of interest
  for (j in 1:length(motifs_of_interest)) {

    # Get current motif
    current_motif <- motifs_of_interest[j]

    # Calculate overlaps for hypergeometric test (upreg)

    # x = upreg and contains motif
    x <- nrow(results[which(
      (results$Label == current_comparison) &
        (results$log2FC > l2fc_max) &
        (results$pvalue < pval_max) &
        (results$adj.pvalue < adjpval_max) &
        (str_detect(string = results$Flankseq, pattern = current_motif))),])


    # m = # motif matches
    m <- nrow(results[which(
      (results$Label == current_comparison) &
        (str_detect(string = results$Flankseq, pattern = current_motif))),])

    # n = # not upreg
    n <- nrow(results[which(results$Label == current_comparison),]) - m

    # k = upreg
    k <- nrow(results[which(
      (results$Label == current_comparison) &
        (results$log2FC > l2fc_max) &
        (results$pvalue < pval_max) &
        (results$adj.pvalue < adjpval_max)),])

    # Calculate hypergeometric pvalue
    pvalue <- dhyper(x, m, n, k)

    # Add to the motif enrichment table
    motif_enrichment_table <- rbind(
      data.frame(
        comparison = current_comparison,
        motif = current_motif,
        x = x,
        m = m,
        n = n,
        k = k,
        pvalue = pvalue,
        direction = "Upreg"
      ) , motif_enrichment_table
    )



    # Calculate overlaps for hypergeometric test (downreg)

    # x = downreg and contains motif
    x <- nrow(results[which(
      (results$Label == current_comparison) &
        (results$log2FC < l2fc_min) &
        (results$pvalue < pval_max) &
        (results$adj.pvalue < adjpval_max) &
        (str_detect(string = results$Flankseq, pattern = current_motif))),])


    # m = # motif matches
    m <- nrow(results[which(
      (results$Label == current_comparison) &
        (str_detect(string = results$Flankseq, pattern = current_motif))),])

    # n = # not upreg
    n <- nrow(results[which(results$Label == current_comparison),]) - m

    # k = upreg
    k <- nrow(results[which(
      (results$Label == current_comparison) &
        (results$log2FC < l2fc_min) &
        (results$pvalue < pval_max) &
        (results$adj.pvalue < adjpval_max)),])

    # Calculate hypergeometric pvalue
    pvalue <- dhyper(x, m, n, k)

    # Add to the motif enrichment table
    motif_enrichment_table <- rbind(
      data.frame(
        comparison = current_comparison,
        motif = current_motif,
        x = x,
        m = m,
        n = n,
        k = k,
        pvalue = pvalue,
        direction = "Downreg"
      ) , motif_enrichment_table
    )


  }



}

# Rename the comparisons
motif_enrichment_table$timepoint <- str_match(string = motif_enrichment_table$comparison, pattern = "PIC_(\\S+)-Control")[,2]



# Adjust for multiple testing
motif_enrichment_table$adj.pvalue <- p.adjust(p = motif_enrichment_table$pvalue, method = "fdr")

# Add logp column
motif_enrichment_table$logp <- -1 * log(motif_enrichment_table$pvalue) / log(10)



# Split table into up and downreg
motif_enrichment_table_up <- motif_enrichment_table[which(motif_enrichment_table$direction == "Upreg"),]

# Filter out RxxS/T from upreg table
motif_enrichment_table_up <- motif_enrichment_table_up[which(!(motif_enrichment_table_up$motif == "....R..........")),]

motif_enrichment_table_down <- motif_enrichment_table[which(motif_enrichment_table$direction == "Downreg"),]

# Filter out downreg  motifs that aren't RxxS/T
motif_enrichment_table_down <- motif_enrichment_table_down[which(motif_enrichment_table_down$motif == "....R.........."),]

# Put them back together again
motif_enrichment_table <- rbind(motif_enrichment_table_up, motif_enrichment_table_down)

# Flip the sign on logp for downreg motifs
motif_enrichment_table[which(motif_enrichment_table$direction == "Downreg"), "logp"] <- -1 * motif_enrichment_table[which(motif_enrichment_table$direction == "Downreg"), "logp"]


# Cast upreg motifs to wide format
motifs_wide <- dcast(data = motif_enrichment_table, formula = motif ~ comparison, value.var = "logp")


# Sort by 2h logp
motifs_wide <- motifs_wide[order(motifs_wide$`PIC_02h-Control`, decreasing = FALSE),]

# Get motif order
motif_order <- motifs_wide$motif



# Plot as tile heatmap
ggplot(data = motif_enrichment_table, aes(x = motif_enrichment_table$timepoint, y = motif_enrichment_table$motif, fill = motif_enrichment_table$logp)) + geom_tile() + scale_y_discrete(limits = motif_order, expand = c(0,0)) + scale_x_discrete(limits = c("05m", "30m", "02h", "24h"), expand = c(0,0)) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "-LogP") + theme_bw() + xlab(label = "Comparison") + ylab(label = "Motif") + ggtitle(label = "Top Enriched Motifs")
ggsave(filename = "20230301-AK10-ph-top-motifs-heatmap.pdf", width= 4, height = 4)







# Plot sequence logos for each motif of interest at 2h
for (i in 1:length(motifs_of_interest)) {

  # Get current motif
  current_motif <- motifs_of_interest[i]

  # Check if it's RxxS
  if (current_motif == "....R..........") {

    # Get flankseqs that match in downreg data
    sequences <- results[which(
      (results$log2FC < l2fc_min) &
        (results$pvalue < pval_max) &
        (results$adj.pvalue < adjpval_max) &
        (results$Label == "PIC_02h-Control") &
        (str_detect(string = results$Flankseq, pattern = current_motif))), "Flankseq"]

    ggseqlogo(data = sequences, method = 'bits')
    ggsave(filename = paste("logo-", i, ".pdf", sep = ""), width = 4, height = 2)

  } else {

    # Get flankseqs that match in downreg data
    sequences <- results[which(
      (results$log2FC > l2fc_max) &
        (results$pvalue < pval_max) &
        (results$adj.pvalue < adjpval_max) &
        (results$Label == "PIC_02h-Control") &
        (str_detect(string = results$Flankseq, pattern = current_motif))), "Flankseq"]

    ggseqlogo(data = sequences, method = 'bits')
    ggsave(filename = paste("logo-", i, ".pdf", sep = ""), width = 4, height = 2)

  }


}
