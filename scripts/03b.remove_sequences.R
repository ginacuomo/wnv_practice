# 03b.remove_sequences.R

# Author: Gina Cuomo-Dannenburg
# Date: 19/01/2026
# 
# Purpose: Remove problematic sequences from second TempEst analysis
# Input: fasta file that is output of 03a
# Output: the above data but with sequences removed that cause issues in TempEst

# library
library(Biostrings)

# read files
seq <- readDNAStringSet("data-generation/03a_france_wnv_filter.fasta")

# sequences to remove based on TempEst analysis of the seq dataset
remove_seq <- c("PV054409_2023-09-15",
                "OZ261140_2022-09-30")
# note: added this info to the metadata notes

# check the sequences are in the original and there's no typos
sum(remove_seq %in% names(seq)) == length(remove_seq)

# remove these sequences
filtered_seq <- seq[!names(seq) %in% remove_seq]
length(seq) - length(filtered_seq) == length(remove_seq)

# save the output
writeXStringSet(filtered_seq, "data-generation/03b_france_wnv_filter.fasta")
