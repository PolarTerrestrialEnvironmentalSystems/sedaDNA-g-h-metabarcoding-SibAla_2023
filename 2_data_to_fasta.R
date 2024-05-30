# TAXALOSS Project - SibAl55-10 matches to GBIF
#  Script from Thomas Boehmer - modified by Jeremy Courtin
# Script last update - 16.03.2023

# Layout of the script:
# 1 - Import modified and manually checked datafile
# 2 - Remove sequences with weird nucleotides
# 3 - Measure coverage
# 4 - Set the columns for fasta file
# 5 - Convert to fasta file

###############################################################################
# GENERAL OPTION - SETTING UP SCRIPT
###############################################################################
rm(list = ls()) # Remove all the objects we created so far.

# Set up options for the rest of the script
options(stringsAsFactors=FALSE) # state that everything is reads as character 

# # Set working directory 


# Load all needed packages
library(readxl)
library(stringr)
library(tidyverse)
library(seqinr)

###############################################################################
# 1 - Import modified and manually checked datafile
###############################################################################
# /!\ New columns with ncbi_species, ncbi_species_id, ncbi_genus, ncbi_genus_id, ncbi_family and ncbi_family_id 
# /!\ manually added to the csv and checked against NCBI for missing information on the ncbi_acc columns
taxatab <- read_delim("SibAl55-10_final_manually_checked_for_fasta.csv", delim = ",", col_names = T) %>% filter(!is.na(DNA_SEQ))
taxatab %>% filter(is.na(ncbi_id_acc))

###############################################################################
# 2 - Remove sequences with weird nucleotides
###############################################################################
taxatab_corrupt <- taxatab[
  (grepl("b", taxatab$DNA_SEQ) | 
     grepl("d", taxatab$DNA_SEQ) |
     grepl("e", taxatab$DNA_SEQ) | 
     grepl("f", taxatab$DNA_SEQ) | 
     grepl("h", taxatab$DNA_SEQ) | 
     grepl("i", taxatab$DNA_SEQ) | 
     grepl("j", taxatab$DNA_SEQ) |
     grepl("k", taxatab$DNA_SEQ) | 
     grepl("l", taxatab$DNA_SEQ) | 
     grepl("m", taxatab$DNA_SEQ) | 
     grepl("n", taxatab$DNA_SEQ) | 
     grepl("o", taxatab$DNA_SEQ) | 
     grepl("p", taxatab$DNA_SEQ) | 
     grepl("q", taxatab$DNA_SEQ) | 
     grepl("r", taxatab$DNA_SEQ) | 
     grepl("s", taxatab$DNA_SEQ) | 
     grepl("u", taxatab$DNA_SEQ) | 
     grepl("v", taxatab$DNA_SEQ) | 
     grepl("w", taxatab$DNA_SEQ) | 
     grepl("x", taxatab$DNA_SEQ) | 
     grepl("y", taxatab$DNA_SEQ) | 
     grepl("z", taxatab$DNA_SEQ)), ]

taxatab_notcorrupt <- taxatab[
  !(grepl("b", taxatab$DNA_SEQ) | 
      grepl("d", taxatab$DNA_SEQ) |
      grepl("e", taxatab$DNA_SEQ) | 
      grepl("f", taxatab$DNA_SEQ) | 
      grepl("h", taxatab$DNA_SEQ) | 
      grepl("i", taxatab$DNA_SEQ) | 
      grepl("j", taxatab$DNA_SEQ) |
      grepl("k", taxatab$DNA_SEQ) | 
      grepl("l", taxatab$DNA_SEQ) | 
      grepl("m", taxatab$DNA_SEQ) | 
      grepl("n", taxatab$DNA_SEQ) | 
      grepl("o", taxatab$DNA_SEQ) | 
      grepl("p", taxatab$DNA_SEQ) | 
      grepl("q", taxatab$DNA_SEQ) | 
      grepl("r", taxatab$DNA_SEQ) | 
      grepl("s", taxatab$DNA_SEQ) | 
      grepl("u", taxatab$DNA_SEQ) | 
      grepl("v", taxatab$DNA_SEQ) | 
      grepl("w", taxatab$DNA_SEQ) | 
      grepl("x", taxatab$DNA_SEQ) | 
      grepl("y", taxatab$DNA_SEQ) | 
      grepl("z", taxatab$DNA_SEQ)), ]

# Check
taxatab[str_detect(taxatab$DNA_SEQ, "y"), ]
taxatab_notcorrupt[str_detect(taxatab_notcorrupt$DNA_SEQ, "y"), ]

# Save as table before fasta file
write_delim(taxatab_notcorrupt, "output/2023_03_14_Table_before_fasta_file_final.csv", delim = ",")

###############################################################################
# 3 - Measure coverage
###############################################################################
taxatab_noseq <- read_delim("SibAl55-10_final_manually_checked_for_fasta.csv", delim = ",", col_names = T) %>% filter(is.na(DNA_SEQ))
coverage <- rbind(taxatab_notcorrupt, taxatab_noseq) %>% mutate(coverage_ext = ifelse(is.na(DNA_SEQ), 0, 1))
coverage %>% select(DNA_SEQ) %>% distinct() %>% filter(!is.na(DNA_SEQ))
coverage %>% select(ncbi_family_acc) %>% distinct() # 233 total families
coverage %>% filter(coverage_ext == 1) %>% select(ncbi_family_acc) %>% distinct() # 223 families in database -> 95.7%
coverage %>% select(ncbi_genus_acc) %>% distinct() # 1059 total genus
coverage %>% filter(coverage_ext == 1) %>% select(ncbi_genus_acc) %>% distinct() # 947 genus in database -> 89.4%
coverage %>% select(ncbi_species_acc) %>% distinct() # 4849 total species
coverage %>% filter(coverage_ext == 1) %>% select(ncbi_species_acc) %>% distinct() # 3398 species in database -> 70.1%
taxatab_notcorrupt %>% select(DNA_SEQ) %>% distinct()
taxatab_notcorrupt %>% select(ncbi_species_acc) %>% distinct()

###############################################################################
# 4 - Set the columns for fasta file
###############################################################################
taxatab2 <- taxatab_notcorrupt %>% mutate(final_name = ncbi_species_acc, species_id = ncbi_species_id, GBIF_tax_status = "accepted") %>% 
  select(DNA_SEQ, final_name, ncbi_species, GBIF_id_acc, ncbi_family_acc, ncbi_family, species_id, ncbi_family_id, ncbi_species_id, GBIF_tax_status, ncbi_genus_acc, ncbi_genus, ncbi_genus_id, db) 
colnames(taxatab2) <- c("sequence", "final_name", "species_name", "GBIF_id_accepted", "family_final","family_name", "taxid", "family", "species", "GBIF_status", "genus_final", "genus_name", "genus", "ori_db")
taxatab_test <- taxatab2 %>% mutate(rank = "species") %>% 
  select(sequence, final_name, family_name, family, family_final, genus_name, genus, genus_final, species_name, species, taxid, rank, ori_db) %>% mutate(ori_db = gsub("[.]", "",ori_db))
taxatab_test[is.na(taxatab_test)] <- "no_entrie"
taxatab_test %>% filter(genus_final == "no_entrie")
taxatab2 <- taxatab_test

###############################################################################
# 5 - Convert to fasta file
###############################################################################
# create table with unique sequences and assigned identifier:
taxatab_uniseq <- unique(taxatab$DNA_SEQ)

# Give you sequence number
sequencetypenumbers <- data.frame(NUC_SEQ_arc = taxatab_uniseq, id=paste0("Sib5510Seq_", sprintf("%05d", seq(1,length(taxatab_uniseq),1))), number_long = sprintf("%04d", seq(1,length(taxatab_uniseq),1)))
write_delim(sequencetypenumbers, "output/SibAl55-10_final_for_fasta_sequencetypes_sequencenumbers_siberia55_occ10_KS.csv", delim = ",")

# check if TaxIDs contain NA:
taxatab2[which(is.na(taxatab2$taxid)), ]

# add Identifier to the datatable
taxatab_final <- data.frame(taxatab2, id=paste0("Sib5510Seq_", sprintf("%05d", seq(1,dim(taxatab2)[1]))))

# extract columns for fasta-annotations:
annot_cols <- taxatab_final[ ,-1]

# merge all annotation columns together:
annot_df <- NULL

for(i in 1:dim(annot_cols)[1]){
  
  print(paste0(i,"/",dim(annot_cols)[1]))
  
  subseq <- annot_cols[i, ]

  subseq_annot_merge <-   str_c(    paste0(names(subseq),"=", subseq)  , collapse = "; ")
  
  annot_df <- c(annot_df, subseq_annot_merge)
  
}

# format sequences as list:
seqs <- as.list(dplyr::pull(taxatab_final, sequence))

# make list of Identifiers:
seqnames <- dplyr::pull(taxatab_final, id)

# bind Identifiers and annotations together:
seqnames_annot <- paste(seqnames, annot_df) 

# write fasta-file - if working with windows need to replace line jumps to LF only before trying to build obitools.
write.fasta(seqs, seqnames_annot, file.out="2023-03-17_Sib5510Seq_sequencetable_siberia55_occ10_clean_nodots_withrank_final.fasta",
            open = "w", as.string = FALSE)


