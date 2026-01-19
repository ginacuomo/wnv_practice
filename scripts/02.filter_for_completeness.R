# 1.edit_seq_names.R

# Author: Gina Cuomo-Dannenburg
# Date: 14/01/2026
# 
# Purpose: Edit the name of the sequences in the fasta file to include the collection date
# Input: fasta file that has been aligned using mafft and manual alignment, and then updated sequence names added
# Output: the above data but with sequences remove that contain less than 80% complete data

# library
library(Biostrings)

seqs <- readDNAStringSet("data-generation/france_wnv_aligned_renamed.fasta")
seq_char <- as.character(seqs)

# split into codons
codons <- lapply(seq_char, function(x) {
  n <- nchar(x)
  starts <- seq(1, n - 2, by = 3)
  substring(x, starts, starts + 2)
})

# define what we want to count as informative codons of information
is_informative <- function(codon) {
  grepl("^[ACGT]{3}$", codon)
}

# TODO: write some tests for this function

# count missing and informative codons
codon_stats <- lapply(codons, function(cds) {
  informative <- sum(sapply(cds, is_informative))
  total <- length(cds)
  missing <- total - informative
  
  data.frame(
    total_codons = total,
    informative_codons = informative,
    missing_codons = missing,
    fraction_missing = missing / total
  )
})

codon_stats <- do.call(rbind, codon_stats)
rownames(codon_stats) <- names(seqs)

# quickly look at distribution of missing data
hist(codon_stats$fraction_missing, breaks = 50)

# set a threshold
threshold <- 0.1 # keep sequences that are at least 20% complete

omit_sequences <- codon_stats %>%
  filter(fraction_missing > threshold) %>%
  rownames()
