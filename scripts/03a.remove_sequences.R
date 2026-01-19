# 03a.remove_sequences.R

# Author: Gina Cuomo-Dannenburg
# Date: 14/01/2026
# 
# Purpose: Remove problematic sequences from TempEst analysis
# Input: fasta file that has been aligned using mafft and manual alignment with filters for completeness.
# Output: the above data but with sequences removed that cause issues in 

# library
library(Biostrings)

# read files
seq <- readDNAStringSet("data-generation/02_france_wnv_aligned_renamed.fasta")

# sequences to remove based on TempEst analysis of the seq dataset
remove_seq <- c("MT863559_2015-10-03",
                "DQ786572_2004-10-18",
                "DQ786573_2004-10-21")
# note: added this info to the metadata notes

# check the sequences are in the original and there's no typos
sum(remove_seq %in% names(seq)) == length(remove_seq)

# remove these sequences
filtered_seq <- seq[!names(seq) %in% remove_seq]
length(seq) - length(filtered_seq) == length(remove_seq)

# save the output
writeXStringSet(filtered_seq, "data-generation/03a_france_wnv_filter.fasta")
