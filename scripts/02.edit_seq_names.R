# 02.edit_seq_names.R

# Author: Gina Cuomo-Dannenburg
# Date: 14/01/2026
# 
# Purpose: Edit the name of the sequences in the fasta file to include the collection date
# Input: fasta file that has undergone mafft alignment and then manually aligning in AliView1.30

# libraries
library(Biostrings)
library(tidyverse)

# read input files
seq <- readDNAStringSet("data-generation/01_france_wnv_threshold.fasta")
meta <- read.csv("data/metadata.csv")

# obtain accession numbers -- all are version .1 so omit this from the name here
seq_accession <- sub("[ |].*$", "", names(seq))
seq_accession <- sub("\\..*$", "", seq_accession)

# prepare metadata
meta <- meta %>%
  filter(Accession %in% seq_accession)
# reorder before we join
meta <- meta[match(seq_accession, meta$Accession), ]
# check that we are not missing things
stopifnot(all(meta$Accession == seq_accession))

# alter collection date to be of the form I want
meta$Collection_Date <- as.Date(meta$Collection_Date, format = "%d/%m/%Y")

# generate new names
new_names <- paste0(meta$Accession, "_", meta$Collection_Date)
any(duplicated(new_names)) # check these are unique

#  rename the FASTA sequences
names(seq) <- new_names

# output the file
writeXStringSet(seq, "data-generation/02_france_wnv_aligned_renamed.fasta")
