# TAXALOSS Project - SibAl55-10 matches to GBIF
# Script from Jeremy Courtin
# Script last update - 16.03.2023

# Layout of the script: 
# 0 - GENERAL OPTION - SETTING UP SCRIPT
# 1 - Download both GBIF and NCBI taxonomies
# 2 - Load GBIF data from the SibAl55-10 area
# 3 - Add NCBI taxonomy to GBIF
# 4 - Load in EMBL, arctborbryo and phylonorway sequence databases
# /!\ manual check missing phylum with NCBI /!\
# 5 - Merge clean Sequence database to GBIF
# /!\ Retry adding NCBI taxonomy to GBIF names! /!\
# /!\ manual check of the list for taxa present in arctic and phylo norway /!\ -> manual step
# 6 - Remove contaminant phase -> from manual .csv
# 7 - Check for synonyms! Keep only accepted names as unique species
# 8 - Final check names directly from the taxdir latest NCBI download
# 9 - Control and save 

###############################################################################
# 0 - GENERAL OPTION - SETTING UP SCRIPT
###############################################################################
rm(list = ls()) # Remove all the objects we created so far.

# Set up options for the rest of the script
options(stringsAsFactors=FALSE) # state that everything is reads as character 

# # Set working directory 


# Load all needed packages
library(tidyverse)
library(taxonbridge)
library(taxize)
library(myTAI)
library(crayon)
library(CHNOSZ)
library(reshape)
library(readr)
library(plyr)
library(taxreturn)

# latest download of the NCBI taxdir
taxdir="NCBI/taxdump"

###############################################################################
# 1 - Download both GBIF and NCBI taxonomies
###############################################################################
# /!\ This needs to be run on linux (could not work on windows for some reason) /!\
#custom_taxonomy <- load_taxonomies(download_gbif(), download_ncbi(taxonkitpath = "/home/jcourtin/bin/taxonkit"))
#save(custom_taxonomy, file = "input/NCBI_to_GBIF_taxonomies.RData")
load(file = "input/NCBI_to_GBIF_taxonomies.RData")

# deduplicate taxonomy
dedup_tax <- dedupe(custom_taxonomy)

###############################################################################
# 2 - Load GBIF data from the SibAl55-10 area
###############################################################################
load("output/steptophyta_combined.RData")

# put it as tibble and change name
gbif_raw <- streptophyta_combined %>% as_tibble()

# quick check of what is in there
gbif_raw_explor <- gbif_raw %>% select(scientificName, taxonKey, acceptedScientificName, acceptedTaxonKey, family, genus, species, taxonRank) %>% distinct()
gbif_raw_explor %>% select(acceptedTaxonKey) %>% distinct() # 10556 accepted taxon key
gbif_raw_explor %>% select(acceptedTaxonKey, taxonRank) %>% distinct() %>% group_by(taxonRank) %>% summarise(nb = length(acceptedTaxonKey)) %>% arrange(nb)
gbif_raw_explor %>% filter(species > 0) %>% select(species) %>% distinct() # 7841 species names
gbif_raw_explor %>% filter(genus > 0) %>% select(genus) %>% distinct() # 1640 genus names
gbif_raw_explor %>% filter(family > 0) %>% select(family) %>% distinct() # 337 family names

# tidy it and get occurrence values
gbif_id <- gbif_raw  %>% filter(!is.na(taxonKey)) %>% dplyr::group_by(acceptedTaxonKey, family, genus, species) %>% dplyr::mutate(occurence = length(acceptedTaxonKey)) %>% ungroup() %>%
  dplyr::select(taxonKey, acceptedTaxonKey, occurence) %>% distinct() %>% arrange(occurence)  %>% filter(occurence > 9)

###############################################################################
# 3 - Add NCBI taxonomy to GBIF
###############################################################################
dedup_tax %>% arrange(taxonID)
gbif <- gbif_id %>% mutate(taxonKey = as.character(taxonKey)) %>% left_join(dedup_tax, by = c("taxonKey" = "taxonID")) %>% distinct()

# We only keep species
gbif_sp <- gbif %>% filter(taxonRank == "species") # 6951 entries

# Tidy info on gbif
names(gbif_sp)
gbif_tidy <- gbif_sp %>% select(occurence, taxonKey, acceptedTaxonKey, taxonRank, taxonomicStatus, family, genericName, specificEpithet, ncbi_id, ncbi_family, ncbi_family, ncbi_genus, ncbi_species) %>% distinct()

# Add accepted species names to table
gbif_acc <- gbif_sp %>% select(acceptedTaxonKey, taxonRank, taxonomicStatus, family, genericName, specificEpithet, ncbi_id, ncbi_family, ncbi_family, ncbi_genus, ncbi_species) %>% filter(taxonomicStatus == "accepted") %>% distinct() 
names(gbif_acc) <- paste("acc", names(gbif_acc), sep = "_")
gbif_tidy_clean <- gbif_tidy %>% left_join(gbif_acc, by = c("acceptedTaxonKey" = "acc_acceptedTaxonKey")) %>% 
  select(occurence, taxonKey, acceptedTaxonKey, taxonRank, acc_taxonRank, taxonomicStatus, family, genericName, specificEpithet, acc_family, acc_genericName, acc_specificEpithet, 
         ncbi_id, acc_ncbi_id, ncbi_family, ncbi_genus, ncbi_species, acc_ncbi_family, acc_ncbi_genus, acc_ncbi_species)

gbif_tidy_clean %>% filter(taxonomicStatus == "accepted") %>% filter(taxonRank == "species") # 4519 species with accepted species
gbif_tidy_clean %>% filter(!taxonomicStatus == "accepted") %>% dplyr::group_by(taxonomicStatus, acc_taxonRank) %>% dplyr::summarise(dat = length(acc_taxonRank)) %>% 
  dplyr::group_by(acc_taxonRank) %>% dplyr::mutate(sumNA = sum(dat)) 
# here some taxa will not have an accepted name  794 without accepted name and 1638 with accepted name -> 2432 taxa with non accepted name
# A total of 6951 taxa: 4519 accepted species and 2432 non-accepted (1638 with accepted link and 794 without)

# CHECKPOINT GBIF INFORMATION
write_delim(gbif_tidy_clean, "CHECKPOINT_GBIF_SibAl55-10area_10occu_species_NCBIlineage.csv", delim = ",")
gbif_tidy_clean <- read_delim("CHECKPOINT_GBIF_SibAl55-10area_10occu_species_NCBIlineage.csv", delim = ",", col_names = T)

###############################################################################
#  4 - Load in EMBL, arctborbryo and phylonorway sequence databases
###############################################################################
# read embl-database and extracted species list from arcborbryo- and phylonorway -> orinigates from Thomas Boehmer conversion from fasta to datatables:
arc_df   <- read_delim("input/DNA_database/2022-05-05_arcbryo_GH_datatable.csv", delim = ",", col_names = T) %>% mutate(taxid = as.character(taxid))
phylo_df <- read_delim("input/DNA_database/2022-05-05_phylonorway_GH_datatable.csv", delim = ",", col_names = T) %>% mutate(taxid = as.character(taxid))
embl_df  <- read_delim("input/DNA_database/2022-05-11_embl143_GH_datatable.csv", delim = ",", col_names = T) %>% mutate(taxid = as.character(taxid))

# merge arctic and phylonorway database
arc <- arc_df %>% select(Seq, taxid) %>% mutate(arc_db = "arc") %>% distinct() # 2289 length
arc %>% select(taxid, arc_db) %>% filter(arc_db == "arc") %>% distinct() # 2099 unique taxid
arc %>% select(Seq) %>% distinct() # 1055 unique sequences
phylo <- phylo_df %>% select(Seq, taxid) %>% mutate(phylo_db = "phylo") %>% distinct() # 1633 length
phylo %>% select(taxid, phylo_db) %>% filter(phylo_db == "phylo") %>% distinct() # 1575 unique taxid
phylo %>% select(Seq) %>% distinct() # 814 unique sequences

arc_phylo <- full_join(arc, phylo) # 3233 length
arc_phylo %>% select(taxid) %>% distinct() # 2920 unique taxid
arc_phylo %>% select(Seq) %>% distinct() # 1381 unique sequences
arc_phylo %>% filter(arc_db == "arc") # 2289 length
arc_phylo %>% select(taxid, arc_db) %>% filter(arc_db == "arc") %>% distinct() # 2099 unique taxid
arc_phylo %>% select(Seq, arc_db) %>% filter(arc_db == "arc") %>% distinct() # 1055 unique sequences
arc_phylo %>% filter(phylo_db == "phylo") # 1633 length
arc_phylo %>% select(taxid, phylo_db) %>% filter(phylo_db == "phylo") %>% distinct() # 1575 unique taxid
arc_phylo %>% select(Seq, phylo_db) %>% filter(phylo_db == "phylo") %>% distinct() # 814 unique sequences

# Get eukaryota from embl
embl <- embl_df %>% filter(superkingdom_name == "Eukaryota") %>% select(Seq, taxid, rank, family, family_name, genus, genus_name, species_name) %>% mutate(embl_db = "embl") %>% distinct() # 90504 length
embl %>% select(taxid, embl_db) %>% filter(embl_db == "embl") %>% distinct() # 75540 unique taxid
embl %>% select(Seq) %>% distinct() # 28193 unique sequences

# Merge all databases
seq_db <- embl %>% full_join(arc_phylo) %>% distinct() %>% mutate(taxid = as.character(taxid)) # 91687
seq_db %>% filter(is.na(embl_db))  %>% distinct() # 1183 entries not in embl but in arc / phylo
seq_db %>% filter(arc_db == "arc") # 2289 length
seq_db %>% filter(phylo_db == "phylo") # 1633 length
seq_db %>% filter(embl_db == "embl") # 90504 length

# correct missing lineage
arc_miss_info <- seq_db %>% filter(arc_db == "arc") %>% filter(is.na(rank)) %>% select(Seq, taxid, embl_db, arc_db, phylo_db) %>% left_join(arc_df) %>% 
  select(Seq, taxid, embl_db, arc_db, phylo_db, rank, family, family_name, genus, genus_name, species, species_name) %>% distinct()
phylo_miss_info <- seq_db %>% filter(phylo_db == "phylo") %>% filter(is.na(rank)) %>% select(Seq, taxid, embl_db, arc_db, phylo_db) %>% left_join(phylo_df) %>% 
  select(Seq, taxid, embl_db, arc_db, phylo_db, rank, family, family_name, genus, genus_name, species, species_name) %>% distinct()
embl_miss_info <- seq_db %>% filter(embl_db == "embl") %>% filter(is.na(rank)) %>% select(Seq, taxid, embl_db, arc_db, phylo_db) %>% left_join(embl_df) %>% 
  select(Seq, taxid, embl_db, arc_db, phylo_db, rank, family, family_name, genus, genus_name, species, species_name) %>% distinct()

miss <- rbind(arc_miss_info, phylo_miss_info, embl_miss_info) %>% select(-embl_db, -arc_db, -phylo_db) %>% distinct()
names(miss) <- paste("corr", names(miss), sep = "_")
miss_clean <- miss %>% mutate(taxid = as.character(corr_taxid), Seq = corr_Seq) %>% select(-corr_taxid) %>% distinct()

# re-add missing taxa to seq_db
seq_db_corrected <- left_join(seq_db, miss_clean) %>% distinct() # keep only species -> 91693
seq_db_to_check <- seq_db_corrected %>% mutate(rank = ifelse(is.na(corr_rank), rank, corr_rank),
                                            family = ifelse(is.na(corr_family), family, corr_family),
                                            family_name = ifelse(is.na(corr_family_name), family_name, corr_family_name),
                                            genus = ifelse(is.na(corr_genus), genus, corr_genus),
                                            genus_name = ifelse(is.na(corr_genus_name), genus_name, corr_genus_name), 
                                            species_name = ifelse(is.na(corr_species_name), species_name, corr_species_name)) %>% 
  select(-corr_rank, -corr_family, -corr_family_name, -corr_genus, -corr_genus_name, -corr_species, -corr_species_name, -corr_Seq) %>% distinct()

# /!\ manual check overlap names /!\ ->need to correct!
check_miss_F <- seq_db_to_check %>% dplyr::group_by(Seq, taxid) %>% dplyr::mutate(occ = n_distinct(family)) %>% arrange(desc(taxid)) %>% filter(occ > 1) # at family level -> 3 duplicates if duplicate, check on NCBI the taxid and keep the name from NCBI
to_remove_F <- check_miss_F %>% filter(family_name == "Chenopodiaceae")
check_miss_G <- seq_db_to_check %>% dplyr::group_by(Seq, taxid) %>%  dplyr::mutate(occ = n_distinct(genus)) %>% arrange(desc(taxid)) %>% filter(occ > 1) # at genus level -> 2 duplicates if duplicate, check on NCBI the taxid and keep the name from NCBI
to_remove_G <- check_miss_G %>% filter(grepl('Botrychium|Hieracium', genus_name))
check_miss_S <- seq_db_to_check %>%  dplyr::group_by(Seq, taxid) %>% dplyr::mutate(occ = n_distinct(species_name)) %>% arrange(desc(taxid)) %>% filter(occ > 1) # at species level -> 1 duplicate if duplicate, check on NCBI the taxid and keep the name from NCBI
to_remove_S <- check_miss_S %>% filter(grepl('Diphasiastrum complanatum', species_name))
to_remove <- rbind(to_remove_F, to_remove_G, to_remove_S) %>% distinct()

seq_db_final <- seq_db_to_check %>% left_join(to_remove) %>% filter(is.na(occ)) %>% select(-occ) %>% distinct() # 91687 -> back to original value!
seq_db_final %>% select(Seq, taxid) %>% distinct() # 91687
seq_db_final %>% select(Seq, taxid, embl_db, arc_db, phylo_db) %>% distinct() # 91687
seq_db_final %>% select(rank) %>% distinct()

# Merge with NCBI taxonomy and only keep Streptophyta
NCBI_taxonomy <- dedup_tax %>% select(ncbi_id, ncbi_phylum) %>% distinct()
seq_db_clean <- left_join(seq_db_final, NCBI_taxonomy, by = c("taxid" = "ncbi_id")) %>% distinct()
seq_db_clean %>% filter(is.na(ncbi_phylum)) # still 829 taxa without phylum
seq_db_clean %>% filter(ncbi_phylum == "Streptophyta") # 86699 Streptophyta
seq_db_clean %>% filter(Seq == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg") # Rosaceae Dryas sequence

###############################################################################
# /!\ manual check missing phylum with NCBI /!\
###############################################################################
tax_counter1 <- seq_db_clean %>% filter(is.na(ncbi_phylum)) %>% select(species_name) %>% distinct() # 723 species name to check ....
tax_counter <- sort(unique(tax_counter1$species_name))

##########################################
taxa_lineage <- data.frame(db_taxon_name=NA, rank=NA, NCBI_TaxID=NA, NCBI_kingdom_name=NA, NCBI_kingdom_id=NA, NCBI_phylum_name=NA, NCBI_phylum_id=NA, NCBI_class_name=NA, NCBI_class_id=NA, 
                           NCBI_order_name=NA, NCBI_order_id=NA, NCBI_family_name=NA, NCBI_family_id=NA, NCBI_genus_name=NA, NCBI_genus_id=NA, NCBI_species_name=NA, NCBI_species_id=NA)

for(i in 1:length(tax_counter)){
  
  cat(red("#################### ",paste0(i,"/",length(tax_counter))," ####################","\n"))
  
  taxa_lineage[i,"db_taxon_name"] <- tax_counter[i]
  
  try(
    suppressWarnings(suppressMessages({
      
      taxdata <- taxonomy(organism=tax_counter[i], db="ncbi", output="classification")                             
      
      
      if(dim(taxdata)[2] > 1){
        
        taxa_lineage[i,"NCBI_TaxID"]        <-  as.numeric(tail(taxdata, n=1)$id)
        taxa_lineage[i,"rank"]              <-  tail(taxdata, n=1)$rank
        
        taxa_lineage[i,"NCBI_kingdom_name"] <- ifelse("kingdom" %in% taxdata$rank, filter(taxdata, rank =="kingdom")$name, "NA")
        taxa_lineage[i,"NCBI_kingdom_id"]   <- ifelse("kingdom" %in% taxdata$rank, filter(taxdata, rank =="kingdom")$id, "NA")
        
        taxa_lineage[i,"NCBI_phylum_name"]  <- ifelse("phylum" %in% taxdata$rank, filter(taxdata, rank =="phylum")$name, "NA")
        taxa_lineage[i,"NCBI_phylum_id"]    <- ifelse("phylum" %in% taxdata$rank, filter(taxdata, rank =="phylum")$id, "NA")
        
        taxa_lineage[i,"NCBI_class_name"]   <- ifelse("class" %in% taxdata$rank, filter(taxdata, rank =="class")$name, "NA")
        taxa_lineage[i,"NCBI_class_id"]     <- ifelse("class" %in% taxdata$rank, filter(taxdata, rank =="class")$id, "NA")
        
        taxa_lineage[i,"NCBI_order_name"]   <- ifelse("order" %in% taxdata$rank, filter(taxdata, rank =="order")$name, "NA")  
        taxa_lineage[i,"NCBI_order_id"]     <- ifelse("order" %in% taxdata$rank, filter(taxdata, rank =="order")$id, "NA")  
        
        taxa_lineage[i,"NCBI_family_name"]  <- ifelse("family" %in% taxdata$rank, filter(taxdata, rank =="family")$name, "NA") 
        taxa_lineage[i,"NCBI_family_id"]    <- ifelse("family" %in% taxdata$rank, filter(taxdata, rank =="family")$id, "NA") 
        
        taxa_lineage[i,"NCBI_genus_name"]   <- ifelse("genus" %in% taxdata$rank, filter(taxdata, rank =="genus")$name, "NA") 
        taxa_lineage[i,"NCBI_genus_id"]     <- ifelse("genus" %in% taxdata$rank, filter(taxdata, rank =="genus")$id, "NA") 
        
        taxa_lineage[i,"NCBI_species_name"] <- ifelse("species" %in% taxdata$rank, filter(taxdata, rank =="species")$name, "NA")
        taxa_lineage[i,"NCBI_species_id"]   <- ifelse("species" %in% taxdata$rank, filter(taxdata, rank =="species")$id, "NA")
        
      }else{
        
        taxa_lineage[i,"NCBI_TaxID"]        <- "NA"
        taxa_lineage[i,"NCBI_TaxID_match"]  <- "NA" 
        taxa_lineage[i,"NCBI_TaxID_name"]   <- "NA"
        taxa_lineage[i,"NCBI_TaxID_ULR"]    <- "NA"
        
        taxa_lineage[i,"NCBI_kingdom_name"] <- "NA"
        taxa_lineage[i,"NCBI_kingdom_id"]   <- "NA"
        
        taxa_lineage[i,"NCBI_phylum_name"]  <- "NA"
        taxa_lineage[i,"NCBI_phylum_id"]    <- "NA"
        
        taxa_lineage[i,"NCBI_class_name"]   <- "NA"
        taxa_lineage[i,"NCBI_class_id"]     <- "NA"
        
        taxa_lineage[i,"NCBI_order_name"]   <- "NA"
        taxa_lineage[i,"NCBI_order_id"]     <- "NA"
        
        taxa_lineage[i,"NCBI_family_name"]  <- "NA"
        taxa_lineage[i,"NCBI_family_id"]    <- "NA"
        
        taxa_lineage[i,"NCBI_genus_name"]   <- "NA"
        taxa_lineage[i,"NCBI_genus_id"]     <- "NA"
        
        taxa_lineage[i,"NCBI_species_name"] <- "NA"
        taxa_lineage[i,"NCBI_species_id"]   <- "NA"
        
      } # end else
      
      # -------------------------------------------------------------------------------------------------
      
    })), silent = TRUE)}

##########################################
# Add phylum information to sequence database
upd_phyl <- taxa_lineage %>% select(NCBI_TaxID, NCBI_species_name, NCBI_phylum_name, NCBI_family_name) %>% distinct() 
upd_phyl_miss_taxa <- tax_counter1 %>% left_join(upd_phyl, by = c("species_name" = "NCBI_species_name")) %>% filter(!is.na(species_name))

################
seq_db_corr <- seq_db_clean %>% left_join(upd_phyl_miss_taxa) %>% mutate(NCBI_fam = ifelse(is.na(family_name), NCBI_family_name, family_name), NCBI_phyl = ifelse(is.na(ncbi_phylum), NCBI_phylum_name, ncbi_phylum))
fam_phyl <- seq_db_corr %>% select(NCBI_fam, NCBI_phyl) %>% distinct()
length(unique(fam_phyl$NCBI_fam)) # 924
corr_phyl <- fam_phyl %>% filter(!is.na(NCBI_phyl)) # 918 -> still 6 families missing but maybe we can manually check them then

seq_db_corr1 <- seq_db_corr %>% select(-NCBI_phyl) %>% left_join(corr_phyl) %>% mutate(family_name = NCBI_fam, ncbi_phylum = NCBI_phyl) %>% select(-NCBI_TaxID, -NCBI_phylum_name, -NCBI_family_name, -NCBI_fam, -NCBI_phyl)
seq_db_corr1 %>% filter(is.na(ncbi_phylum)) %>% select(family_name) # Only Perssoniellaceae is Streptophyta!
seq_db_clean_corr <- seq_db_corr1 %>% mutate(ncbi_phylum = ifelse(family_name == "Perssoniellaceae", "Streptophyta", ncbi_phylum)) %>% distinct() # 84595 again

seq_db_clean %>% filter(ncbi_phylum == "Streptophyta") # before 86699 
seq_db_clean_corr %>% filter(ncbi_phylum == "Streptophyta") # now 87435

seq_db_clean_plant <- seq_db_clean_corr %>% mutate(embl_db = ifelse(is.na(embl_db), "", embl_db), arc_db = ifelse(is.na(arc_db), "", arc_db), phylo_db = ifelse(is.na(phylo_db), "", phylo_db)) %>%
  mutate(db = paste(embl_db, arc_db, phylo_db, sep = ".")) %>% filter(ncbi_phylum == "Streptophyta") %>% select(-ncbi_phylum , -embl_db, -arc_db, -phylo_db) %>% distinct() # 87435 Steptophyta

# get only sequence with info for species
seq_db_clean_sp <- seq_db_clean_plant %>% filter(rank == "species") # 81690
seq_db_clean_sp %>% group_by(db) %>% summarise(length(db))
seq_db_clean_sp %>% filter(grepl("arc", db)) # 2136 length -> 2289 length before
seq_db_clean_sp %>% filter(grepl("phylo", db)) # 1527 length -> 1633 length before
seq_db_clean_sp %>% filter(grepl("embl", db)) # 80701 length -> 90504 length before

seq_db_clean_sp %>% select(species_name) %>% distinct()

#CHECK POINT SEQUENCE DATABASE
write_delim(seq_db_clean_sp, "CHECKPOINT_NCBI_database_clean_species.csv", delim = ",", col_names = T)
seq_db_clean_sp <- read_delim("CHECKPOINT_NCBI_database_clean_species.csv", delim = ",", col_names = T)
seq_db_clean_sp %>% filter(Seq == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg")

###############################################################################
# 5 - Merge clean Sequence database to GBIF
###############################################################################
gbif_tidy_clean <- read_delim("CHECKPOINT_GBIF_SibAl55-10area_10occu_species_NCBIlineage.csv", delim = ",", col_names = T) %>% mutate(taxonKey = as.character(taxonKey))

# re-add GBIF taxonomy to SibAl database
seq_GBIF <- seq_db_clean_sp %>% mutate(family = as.character(family), genus = as.character(genus), taxid = as.character(taxid)) %>% mutate(ncbi_id = taxid) %>% left_join(dedup_tax, by = "ncbi_id") %>%
  select(Seq, ncbi_id, rank, family.x, family_name, genus, genus_name, species_name, db, taxonID) %>% distinct() # 81690

seq_GBIF %>% select(ncbi_id, taxonID) %>% distinct() %>% filter(is.na(taxonID)) # -> 5235 ncbi_id without GBIF entries 
seq_GBIF %>% filter(grepl("arc", db)) %>% select(ncbi_id, taxonID) %>% distinct() %>% filter(is.na(taxonID)) # 291 missing arc entries
seq_GBIF %>% filter(grepl("phylo", db)) %>% select(ncbi_id, taxonID) %>% distinct() %>% filter(is.na(taxonID)) # 11 missing phylo entries

# Add database to GBIF
SibAl <- gbif_tidy_clean %>% left_join(seq_GBIF, by = c("taxonKey" = "taxonID")) %>% distinct() %>% filter(!is.na(occurence)) %>% 
  select(occurence, taxonKey, acceptedTaxonKey, taxonomicStatus, Seq, ncbi_id.x, rank, family.x, family_name, genus, genus_name, species_name, db) %>%
  mutate(family = family.x, ncbi_id = ncbi_id.x) %>% select(-family.x, -ncbi_id.x)

###############################################################################
# /!\ Retry adding NCBI taxonomy to GBIF names! /!
###############################################################################
gbif_raw_names <- gbif_raw %>% filter(taxonomicStatus == "ACCEPTED") %>% select(acceptedTaxonKey, genericName, specificEpithet) %>% distinct() %>% mutate(species_name = paste(genericName, specificEpithet, sep = " "))
GBIF_names <- SibAl %>% filter(is.na(ncbi_id)) %>% select(acceptedTaxonKey) %>% distinct() %>% left_join(select(gbif_raw_names, acceptedTaxonKey, species_name)) %>% distinct() %>% filter(!is.na(species_name))
tax_counter <- sort(unique(GBIF_names$species_name)) # check accepted

###########################
taxa_lineage <- data.frame(db_taxon_name=NA, rank=NA, NCBI_TaxID=NA, NCBI_kingdom_name=NA, NCBI_kingdom_id=NA, NCBI_phylum_name=NA, NCBI_phylum_id=NA, NCBI_class_name=NA, NCBI_class_id=NA, 
                           NCBI_order_name=NA, NCBI_order_id=NA, NCBI_family_name=NA, NCBI_family_id=NA, NCBI_genus_name=NA, NCBI_genus_id=NA, NCBI_species_name=NA, NCBI_species_id=NA)

for(i in 1:length(tax_counter)){
  
  cat(red("#################### ",paste0(i,"/",length(tax_counter))," ####################","\n"))
  
  
  taxa_lineage[i,"db_taxon_name"] <- tax_counter[i]
  
  try(
    suppressWarnings(suppressMessages({
      
      taxdata <- taxonomy(organism=tax_counter[i], db="ncbi", output="classification")                             
      
      
      if(dim(taxdata)[2] > 1){
        
        taxa_lineage[i,"NCBI_TaxID"]        <-  as.numeric(tail(taxdata, n=1)$id)
        taxa_lineage[i,"rank"]              <-  tail(taxdata, n=1)$rank
        
        taxa_lineage[i,"NCBI_kingdom_name"] <- ifelse("kingdom" %in% taxdata$rank, filter(taxdata, rank =="kingdom")$name, "NA")
        taxa_lineage[i,"NCBI_kingdom_id"]   <- ifelse("kingdom" %in% taxdata$rank, filter(taxdata, rank =="kingdom")$id, "NA")
        
        taxa_lineage[i,"NCBI_phylum_name"]  <- ifelse("phylum" %in% taxdata$rank, filter(taxdata, rank =="phylum")$name, "NA")
        taxa_lineage[i,"NCBI_phylum_id"]    <- ifelse("phylum" %in% taxdata$rank, filter(taxdata, rank =="phylum")$id, "NA")
        
        taxa_lineage[i,"NCBI_class_name"]   <- ifelse("class" %in% taxdata$rank, filter(taxdata, rank =="class")$name, "NA")
        taxa_lineage[i,"NCBI_class_id"]     <- ifelse("class" %in% taxdata$rank, filter(taxdata, rank =="class")$id, "NA")
        
        taxa_lineage[i,"NCBI_order_name"]   <- ifelse("order" %in% taxdata$rank, filter(taxdata, rank =="order")$name, "NA")  
        taxa_lineage[i,"NCBI_order_id"]     <- ifelse("order" %in% taxdata$rank, filter(taxdata, rank =="order")$id, "NA")  
        
        taxa_lineage[i,"NCBI_family_name"]  <- ifelse("family" %in% taxdata$rank, filter(taxdata, rank =="family")$name, "NA") 
        taxa_lineage[i,"NCBI_family_id"]    <- ifelse("family" %in% taxdata$rank, filter(taxdata, rank =="family")$id, "NA") 
        
        taxa_lineage[i,"NCBI_genus_name"]   <- ifelse("genus" %in% taxdata$rank, filter(taxdata, rank =="genus")$name, "NA") 
        taxa_lineage[i,"NCBI_genus_id"]     <- ifelse("genus" %in% taxdata$rank, filter(taxdata, rank =="genus")$id, "NA") 
        
        taxa_lineage[i,"NCBI_species_name"] <- ifelse("species" %in% taxdata$rank, filter(taxdata, rank =="species")$name, "NA")
        taxa_lineage[i,"NCBI_species_id"]   <- ifelse("species" %in% taxdata$rank, filter(taxdata, rank =="species")$id, "NA")
        
      }else{
        
        taxa_lineage[i,"NCBI_TaxID"]        <- "NA"
        taxa_lineage[i,"NCBI_TaxID_match"]  <- "NA" 
        taxa_lineage[i,"NCBI_TaxID_name"]   <- "NA"
        taxa_lineage[i,"NCBI_TaxID_ULR"]    <- "NA"
        
        taxa_lineage[i,"NCBI_kingdom_name"] <- "NA"
        taxa_lineage[i,"NCBI_kingdom_id"]   <- "NA"
        
        taxa_lineage[i,"NCBI_phylum_name"]  <- "NA"
        taxa_lineage[i,"NCBI_phylum_id"]    <- "NA"
        
        taxa_lineage[i,"NCBI_class_name"]   <- "NA"
        taxa_lineage[i,"NCBI_class_id"]     <- "NA"
        
        taxa_lineage[i,"NCBI_order_name"]   <- "NA"
        taxa_lineage[i,"NCBI_order_id"]     <- "NA"
        
        taxa_lineage[i,"NCBI_family_name"]  <- "NA"
        taxa_lineage[i,"NCBI_family_id"]    <- "NA"
        
        taxa_lineage[i,"NCBI_genus_name"]   <- "NA"
        taxa_lineage[i,"NCBI_genus_id"]     <- "NA"
        
        taxa_lineage[i,"NCBI_species_name"] <- "NA"
        taxa_lineage[i,"NCBI_species_id"]   <- "NA"
        
      } # end else
      
      # -------------------------------------------------------------------------------------------------
      
    })), silent = TRUE)}

###########################
# Add phylum information to sequence database
updt_GBIF_names <- taxa_lineage %>% filter(rank == "species") %>% select(NCBI_TaxID, NCBI_species_name, NCBI_genus_name, NCBI_family_name) %>% distinct() 
updt_GBIF_names_with_miss <- GBIF_names %>% left_join(updt_GBIF_names, by = c("species_name" = "NCBI_species_name")) %>% filter(!is.na(species_name)) %>% mutate(NCBI_species_name = species_name) %>% select(-species_name)
updt_GBIF_names_with_miss %>% select(NCBI_species_name, acceptedTaxonKey) %>% distinct()

# 
SibAl_corr <- SibAl %>% left_join(updt_GBIF_names_with_miss, by = "acceptedTaxonKey") %>% mutate(GBIF_acc_species = NCBI_species_name) %>% select(-NCBI_species_name) %>% 
  mutate(family_name = ifelse(is.na(family_name), NCBI_family_name, family_name), genus_name = ifelse(is.na(genus_name), NCBI_genus_name, genus_name), species_name = ifelse(is.na(species_name), GBIF_acc_species, species_name)) %>%
  select(-NCBI_family_name, -NCBI_genus_name, -NCBI_TaxID)

write_delim(SibAl_corr, "output/SibAl_database_all.csv", delim = ",")
SibAl_corr <- read_delim("output/SibAl_database_all.csv", delim = ",", col_names = T)
SibAl_corr %>% filter(Seq == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg")

# Re-add extra taxa from arctic and phylonorway
arc_phylo_gbif_pres <- SibAl_corr %>% filter(grepl("arc", db) | grepl("phylo", db)) %>% select(Seq, ncbi_id) %>% distinct() %>% mutate(gbif_pres = "gbif") %>% mutate(ncbi_id = as.character(ncbi_id)) # 2067 NCBI id and sequence type combination present in gbif

seq_GBIF_arc_phylo <- seq_GBIF %>% select(-taxonID) %>% distinct() %>% filter(grepl("arc", db) | grepl("phylo", db)) %>% left_join(arc_phylo_gbif_pres) %>% distinct()
seq_GBIF_arc_phylo %>% select(Seq, ncbi_id) %>% distinct() # 2984 Seq and NCBI id combination
seq_GBIF_arc_phylo %>% select(ncbi_id) %>% distinct() # 2690 unique NCBI_id
seq_GBIF_arc_phylo %>% select(Seq, ncbi_id, gbif_pres) %>% distinct() %>% dplyr::group_by(gbif_pres) %>% dplyr::summarise(length(gbif_pres)) # 2067 seq and NCBI_id present in gbif and 917 absent
seq_GBIF_arc_phylo %>% select(ncbi_id, gbif_pres) %>% distinct() %>% dplyr::group_by(gbif_pres) %>% dplyr::summarise(length(gbif_pres)) # 1827 NCBI_id present and 863 absent

SibAl_extra <- seq_GBIF_arc_phylo %>% filter(is.na(gbif_pres)) # 917 lines
SibAl_extra %>% select(ncbi_id) %>% distinct() # 863 missing NCBI_id
SibAl_extra %>% filter(Seq == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg")
save_phylo_arc <- SibAl_extra %>% select(ncbi_id, db, family_name, genus_name, species_name) %>% distinct()
write_delim(save_phylo_arc, "input/GBIF_absent_phylo_arc_final.csv", delim = ",")

######################################################################
######################################################################
###############################################################################
# /!\ manual check of the list for taxa present in arctic and phylo norway /!\
###############################################################################
######################################################################
######################################################################
# From the saved file: GBIF_absent_phylo_arc_final.csv : check all sequence from arctic and phylonorway database. Are they really absent from the area? Use powo to manually check each entry: 
# add a column to GBIF_absent_phylo_arc_final.csv with keep as a header and put "yes" if we need to keep this sequence or "no" if we do not need it.
# Save the check file under GBIF_absent_phylo_arc_manual_checked.csv and re-import it.
save_phylo_arc_clean <- read_delim("input/GBIF_absent_phylo_arc_manual_checked.csv", delim = ",", col_names = T) %>% filter(keep == "yes") %>% select(db, family_name, genus_name, species_name)

# Add the infomration to SibAl_corr or SibAl_database_all.csv
SibAl_with_extra_arc_phylo <- SibAl_corr %>% full_join(save_phylo_arc_clean)

###############################################################################
# 6 - Remove contaminant phase
###############################################################################
# Import contaminant list - check from Kathleen Stoof-Leichsenring
conta <- read_delim("input/KS_conta_list.csv", delim = ",")
conta_list <- conta %>% filter(contaminants_confirmed == "yes") %>% select(NCBI_taxID) %>% distinct() %>% filter(!is.na(NCBI_taxID))

contaminant_list <- conta_list %>% mutate(NCBI_taxID = as.character((NCBI_taxID))) %>% left_join(custom_taxonomy, by = c("NCBI_taxID" = "ncbi_id"))

# Remove contaminant from list using NCBI id
SibAl_with_extra_arc_phylo_clean <- SibAl_with_extra_arc_phylo %>% filter(!ncbi_id %in% conta_list$NCBI_taxID) %>% 
  setNames(c("GBIF_occurence", "GBIF_id", "GBIF_id_acc", "GBIF_tax_status", "DNA_SEQ", "ncbi_rank", "ncbi_family", "ncbi_genus_id", "ncbi_genus", "ncbi_species", "db", "ncbi_family_id", "ncbi_id", "GBIF_acc_species")) %>%
  select(GBIF_id, GBIF_tax_status, GBIF_id_acc, GBIF_acc_species, GBIF_occurence, ncbi_rank, ncbi_family_id, ncbi_family, ncbi_genus_id, ncbi_genus, ncbi_id, ncbi_species, db, DNA_SEQ) 

SibAl_with_extra_arc_phylo_clean %>% filter(!is.na(DNA_SEQ)) %>% filter(is.na(ncbi_id))

# CHECKPOINT without synonym checks!
write_delim(SibAl_with_extra_arc_phylo_clean, "CHECKPOINT_merge_DNAseq_GBIF_with_synonyms.csv", delim = ",")
SibAl_with_extra_arc_phylo_clean %>% filter(DNA_SEQ == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg")

###############################################################################
# 7 - Check for synonyms! Keep only accepted names as unique species
###############################################################################
# Keep only taxa with accepted GBIF status
accpted <- SibAl_with_extra_arc_phylo_clean %>% filter(GBIF_tax_status == "accepted") %>% select(GBIF_id_acc, ncbi_family_id, ncbi_family, ncbi_genus_id, ncbi_genus, ncbi_id, ncbi_species) %>% 
  setNames(c("GBIF_id_acc", "ncbi_family_id_acc", "ncbi_family_acc", "ncbi_genus_id_acc", "ncbi_genus_acc", "ncbi_id_acc", "ncbi_species_acc")) 

# Arrange and fill missing information for similar accepted taxa names
SibAl_final_names <- SibAl_with_extra_arc_phylo_clean %>% select(-GBIF_acc_species) %>% distinct() %>% left_join(accpted, by = "GBIF_id_acc") %>% distinct() %>% 
  arrange(desc(DNA_SEQ)) %>% dplyr::group_by(GBIF_occurence, ncbi_species_acc, GBIF_id_acc) %>% fill(DNA_SEQ) %>% ungroup() %>%
  arrange(desc(db)) %>% dplyr::group_by(GBIF_occurence, ncbi_species_acc, GBIF_id_acc) %>% fill(db) %>% ungroup() %>%
  mutate(coverage = ifelse(is.na(DNA_SEQ), 0, 1)) %>%
  mutate(GBIF_tax_status = ifelse(is.na(GBIF_id_acc), "accepted", GBIF_tax_status))
SibAl_final_names %>% filter(DNA_SEQ == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg")

# check
SibAl_final_names %>% filter(is.na(GBIF_id_acc))
SibAl_final_names %>% filter(coverage == 1) %>% filter(is.na(ncbi_species_acc))

# Re-add non-accepted names when missing information and only keep accepted taxa
merg_seq <- SibAl_extra %>% select(db, Seq, family.x, family_name, genus, genus_name, ncbi_id, species_name) %>% setNames(c("db", "DNA_SEQ1", "ncbi_family_id_acc1", "ncbi_family", "ncbi_genus_id_acc1", "ncbi_genus", "ncbi_id_acc1", "ncbi_species"))
SibAl_final_names_s <- SibAl_final_names %>% left_join(merg_seq)
SibAl_final_names_t <- SibAl_final_names_s %>% 
  mutate(ncbi_species_acc = ifelse(is.na(GBIF_id_acc), ifelse(GBIF_tax_status == "accepted", ncbi_species, ncbi_species_acc), ncbi_species_acc),
         ncbi_genus_acc = ifelse(is.na(GBIF_id_acc), ifelse(GBIF_tax_status == "accepted", ncbi_genus, ncbi_genus_acc), ncbi_genus_acc),
         ncbi_family_acc = ifelse(is.na(GBIF_id_acc), ifelse(GBIF_tax_status == "accepted", ncbi_family, ncbi_family_acc), ncbi_family_acc),
         ncbi_id_acc = ifelse(is.na(GBIF_id_acc), ifelse(GBIF_tax_status == "accepted", ncbi_id_acc1, ncbi_id_acc), ncbi_id_acc),
         DNA_SEQ = ifelse(is.na(GBIF_id_acc), ifelse(GBIF_tax_status == "accepted", DNA_SEQ1, DNA_SEQ), DNA_SEQ),
         ncbi_genus_id_acc = ifelse(is.na(GBIF_id_acc), ifelse(GBIF_tax_status == "accepted", ncbi_genus_id_acc1, ncbi_genus_id_acc), ncbi_genus_id_acc),
         ncbi_family_id_acc = ifelse(is.na(GBIF_id_acc), ifelse(GBIF_tax_status == "accepted", ncbi_family_id_acc1, ncbi_family_id_acc), ncbi_family_id_acc)) %>% 
  mutate(ncbi_species_acc = ifelse(is.na(ncbi_species_acc), ifelse(coverage == 1, ncbi_species, ncbi_species_acc), ncbi_species_acc),
         ncbi_genus_acc = ifelse(is.na(ncbi_genus_acc), ifelse(coverage == 1, ncbi_genus, ncbi_genus_acc), ncbi_genus_acc),
         ncbi_family_acc = ifelse(is.na(ncbi_family_acc), ifelse(coverage == 1, ncbi_family, ncbi_family_acc), ncbi_family_acc),
         ncbi_id_acc = ifelse(is.na(ncbi_id_acc), ifelse(coverage == 1, ncbi_id_acc1, ncbi_id_acc), ncbi_id_acc),
         ncbi_genus_id_acc = ifelse(is.na(ncbi_genus_id_acc), ifelse(coverage == 1, ncbi_genus_id_acc1, ncbi_genus_id_acc), ncbi_genus_id_acc),
         ncbi_family_id_acc = ifelse(is.na(ncbi_family_id_acc), ifelse(coverage == 1, ncbi_family_id_acc1, ncbi_family_id_acc), ncbi_family_id_acc)) %>%
  select(-DNA_SEQ1, -ncbi_family_id_acc1, -ncbi_family, -ncbi_genus_id_acc1, -ncbi_genus, -ncbi_id_acc1, -ncbi_species, -ncbi_rank, -ncbi_family_id, -ncbi_genus_id) %>% distinct() %>%
  dplyr::group_by(DNA_SEQ, GBIF_occurence, GBIF_id_acc) %>% dplyr::mutate(weird = n_distinct(db)) %>% ungroup() %>% arrange(ncbi_species_acc) %>% filter(GBIF_tax_status == "accepted") %>% select(-coverage, -weird) %>% distinct()

SibAl_final_names_t %>% filter(DNA_SEQ == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg")

###############################################################################
# 8 - Final check names directly from the taxdir latest NCBI download
###############################################################################
# taxon name (eg.scientific name, blast name..., Not nodes!)
ncbi_names=getnames(taxdir)

# read rankedlineage.dmp. 
# The ranks are species, genus, family, order, class, phylum, kingdom, and superkingdom.
rankedlineage=as.data.frame(read_tsv(paste0(taxdir, "/rankedlineage.dmp"),
                                     col_names = c("taxid", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"),
                                     col_types = ("i-c-c-c-c-c-c-c-c-c-")))
# read fullnamelineage.dmp
fullnamelineage=as.data.frame(read_tsv(paste0(taxdir, "/fullnamelineage.dmp"),
                                       col_names = c("taxid", "tax_name", "fullnamelineage"),
                                       col_types = ("i-c-c-")))

ncbi_df=merge(fullnamelineage[c("taxid", "fullnamelineage")], rankedlineage, by = "taxid")
ncbi_df <- ncbi_df %>% mutate(taxid = as.character(taxid))

# attach lineage information to latest the SibAl55-10 data
classified_ncbi <- left_join(SibAl_final_names_t, ncbi_df, by = c("ncbi_id_acc" = "taxid"))
length(unique(classified_ncbi$ncbi_id_acc))

names(classified_ncbi)
classified_ncbi %>% select(tax_name)

# Change the names again
SibAl_final_names_y <- classified_ncbi %>% 
  mutate(ncbi_family_acc = ifelse(is.na(ncbi_family_acc), family, ncbi_family_acc),
         ncbi_genus_acc = ifelse(is.na(ncbi_genus_acc), genus, ncbi_genus_acc),
         ncbi_species_acc = ifelse(is.na(ncbi_species_acc), tax_name, ncbi_species_acc)) %>%
  select(-fullnamelineage, -tax_name, -species, -genus, -family, -family, -order, -class, -phylum, -kingdom, -superkingdom)
SibAl_final_names_y %>% filter(DNA_SEQ == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg")

###############################################################################
# 9 - Control and save 
###############################################################################
# Control if i did not add contaminants again...
conta_sp <- conta %>% filter(contaminants_confirmed == "yes") %>% select(species_name) %>% distinct() %>% filter(!is.na(species_name))
SibAl_final_names_y %>% filter(ncbi_species_acc %in% all_of(conta_sp$species_name))

SibAl_final_names_z <- SibAl_final_names_y %>% filter(!ncbi_species_acc %in% all_of(conta_sp$species_name))

# + check if all added species from arc and phylo are here:
length(unique(save_phylo_arc_clean$species_name)) # 856 sp3ecies to add
SibAl_final_names_z %>% subset(ncbi_species_acc %in% all_of(save_phylo_arc_clean$species_name)) %>% select(ncbi_species_acc) %>% distinct() # 839 added (-contaminants)

#SAVE the SibAl_GBIF_data and hand modify the duplicate accepted ID that have different names!
SibAl_final_names_clean <- SibAl_final_names_z %>% 
  arrange(desc(DNA_SEQ)) %>% group_by(GBIF_occurence, ncbi_species_acc, GBIF_id_acc) %>% fill(DNA_SEQ) %>% ungroup() %>%
  arrange(desc(ncbi_family_id_acc)) %>% group_by(ncbi_family_acc) %>% fill(ncbi_family_id_acc) %>% ungroup() %>%
  arrange(desc(ncbi_genus_id_acc)) %>% group_by(ncbi_family_acc, ncbi_genus_acc) %>% fill(ncbi_genus_id_acc) %>% ungroup() %>%
  arrange(desc(ncbi_id_acc)) %>% group_by(ncbi_family_acc, ncbi_genus_acc, ncbi_species_acc) %>% fill(ncbi_id_acc) %>% ungroup()

SibAl_final_names_clean %>% filter(DNA_SEQ == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg")

ori_NCBI <- SibAl_corr %>% mutate(DNA_SEQ = Seq, GBIF_id_acc = acceptedTaxonKey, GBIF_occurence = occurence, GBIF_id = taxonKey) %>% 
  select(DNA_SEQ, GBIF_id_acc, family, family_name, genus, genus_name, ncbi_id, species_name, db) %>% filter(!is.na(DNA_SEQ)) %>% filter(!is.na(GBIF_id_acc))
SibAl_final_names_clean_final <- SibAl_final_names_clean %>% left_join(ori_NCBI) %>% distinct() %>% 
  mutate(family = ifelse(is.na(ncbi_family_id_acc), family, ncbi_family_id_acc), family_name = ifelse(is.na(ncbi_family_acc), family_name, ncbi_family_acc),
         genus = ifelse(is.na(ncbi_genus_id_acc), genus, ncbi_genus_id_acc), genus_name = ifelse(is.na(ncbi_genus_acc), genus_name, ncbi_genus_acc),
         ncbi_id = ifelse(is.na(ncbi_id_acc), ncbi_id, ncbi_id_acc), species_name = ifelse(is.na(ncbi_species_acc), species_name, ncbi_species_acc)) %>% distinct()
SibAl_final_names_clean_final %>% filter(is.na(ncbi_id_acc)) %>%  filter(!is.na(DNA_SEQ))
names(SibAl_final_names_clean_final)
SibAl_final_names_clean_final %>% filter(DNA_SEQ == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg")

SibAl_final_names_save <- SibAl_final_names_clean_final %>%
  setNames(c("GBIF_id", "GBIF_tax_status", "GBIF_id_acc", "GBIF_occurence", "ncbi_species_id", "db", "DNA_SEQ", "ncbi_family_id_acc", "ncbi_family_acc", "ncbi_genus_id_acc", "ncbi_genus_acc", "ncbi_id_acc", "ncbi_species_acc",
             "ncbi_family_id", "ncbi_family", "ncbi_genus_id", "ncbi_genus", "ncbi_species")) %>% 
  select(DNA_SEQ, GBIF_occurence, GBIF_id_acc, db, ncbi_family_id_acc, ncbi_family_acc, ncbi_genus_id_acc, ncbi_genus_acc, ncbi_id_acc, ncbi_species_acc, ncbi_family_id, ncbi_family, ncbi_genus_id, ncbi_genus, ncbi_species_id, ncbi_species) %>% distinct() %>%
  dplyr::group_by(DNA_SEQ, GBIF_occurence, GBIF_id_acc) %>% dplyr::mutate(duplicates_to_fix = n_distinct(ncbi_species)) %>% ungroup() %>%
  dplyr::group_by(DNA_SEQ, ncbi_species) %>% dplyr::mutate(duplicates_to_fix2 = n()) %>% ungroup()
  
SibAl_final_names_save %>% filter(DNA_SEQ == "atcctgttttatgaaaacaaacaagggtttcagaaagcgcgaataaaaagg")

SibAl_final_names_save_noseq <- select(SibAl_final_names_save, -DNA_SEQ, -db) %>% distinct()

#SAVE the SibAl_GBIF_data
write_delim(SibAl_final_names_clean_final, "output/GBIF_SibAl55-10_database_final.csv", delim = ",")
write_delim(SibAl_final_names_save, "output/GBIF_SibAl55-10_database_final_clean.csv", delim = ",")
write_delim(SibAl_final_names_save_noseq, "output/GBIF_SibAl55-10_taxalist_noseq_final_clean.csv", delim = ",")

# FROM THERE, I check manually and prepare for the fasta -> need to manual check missing family ids, genus ids and species ids as well as the duplicates row! Please check on NCBI and GBIF the best names to keep. 
# Create 2 new colums: synonym_id and add the ncbiIDs of the synonyms of kept taxa and db_mod, add the db where the synonyms are present
# Create new columns with 
write_delim(SibAl_final_names_save, "output/SibAl55-10_final_for_manual_check_for_fasta.csv", delim = ",")
