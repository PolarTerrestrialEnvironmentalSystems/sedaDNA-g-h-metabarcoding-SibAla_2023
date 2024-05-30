# TAXALOSS Project - SibAl55-10 matches to GBIF
# Script from Thomas Boehmer - modified by Jeremy Courtin
# Script last update - 16.03.2023


###############################################################################
# GENERAL OPTION - SETTING UP SCRIPT
###############################################################################
rm(list = ls()) # Remove all the objects we created so far.

# Set up options for the rest of the script
options(stringsAsFactors=FALSE) # state that everything is reads as character 

# # Set working directory 
library(tidyverse)

# read namelist taxatable:
sib55 <- read.csv("input/2022-07-08_sequencetable_taxanames_lineage_siberia55_occ10_final.csv", header=TRUE, sep="\t", check.names="FALSE")

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### Taxa occurrences for "Siberia55"-Region (55-90°N, 50-150°E):

{
  
  # read occurrence table:
  gbif_sib55 <- read.table("input/gbif_taxa_occurences_sib55_region/occurrence.txt", 
                           sep="\t", header=TRUE, comment.char="#", na.strings=".", stringsAsFactors=FALSE, quote="", fill=TRUE)
  
  
  # select only those columns in the dataframe with necessary information:
  keep <- c("gbifID", "occurrenceID", "recordNumber", "individualCount", "organismQuantity", "organismQuantityType", "georeferenceVerificationStatus", "occurrenceStatus", "associatedReferences", "associatedTaxa", "continent", "countryCode", "stateProvince", "county", "municipality", "locality", 
            "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", "verbatimCoordinateSystem", "taxonID", "scientificNameID", "acceptedNameUsageID", "scientificName", "higherClassification", "kingdom", "phylum", "class", "order", "family", "genus", "genericName", "specificEpithet", 
            "taxonRank", "taxonomicStatus", "hasCoordinate", "taxonKey", "acceptedTaxonKey", "kingdomKey", "phylumKey", "classKey", "orderKey", "familyKey", "genusKey", "speciesKey", "species", "acceptedScientificName", "verbatimScientificName")
  
  
  # filter the gbif_sib55 dataframe for only those columns that should be kept:
  gbif_sib55_keep <- gbif_sib55[ ,which(names(gbif_sib55) %in% keep)]
  
  
  # filter taxa in the gbif_sib55_keep dataframe for Streptophyta:
  streptophyta_sib55 <- subset(gbif_sib55_keep, phylum == "Tracheophyta" | phylum == "Bryophyta" | 
                                 phylum == "Marchantiophyta" | phylum == "Anthocerotophyta")
  
  # -------------------------------------------------------------------------------------------------
  
  # calculate number of occurrences for each taxon: 
  
  sib55_IDs <- sort(unique(streptophyta_sib55$scientificName))
  
  gbif_sib55_occurrences <- data.frame(acceptedScientificName=NA, scientificName=NA, family=NA, genus=NA, species=NA, occurrence=NA)
  
  
  for(i in 1:length(sib55_IDs)){
    
    print(paste0(i,"/",length(sib55_IDs)))
    
    taxasubset <- streptophyta_sib55[streptophyta_sib55$scientificName == sib55_IDs[i], ]
    
    accScienNam <- unique(taxasubset$acceptedScientificName)
    
    gbif_sib55_occurrences[i,"acceptedScientificName"] <- accScienNam[which(accScienNam != "")]
    gbif_sib55_occurrences[i,"scientificName"]         <- unique(taxasubset$scientificName)
    gbif_sib55_occurrences[i,"family"]                 <- unique(taxasubset$family) 
    gbif_sib55_occurrences[i,"genus"]                  <- unique(taxasubset$genus)
    gbif_sib55_occurrences[i,"species"]                <- unique(taxasubset$species)
    gbif_sib55_occurrences[i,"occurrence"]             <- dim(taxasubset)[1]
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # filter occurrence dataframe for all taxa with occurrences >= 10:
  gbif_sib55_occurrences10 <- subset(gbif_sib55_occurrences, occurrence >= 10)
  
  # filter streptophyta_sib55 dataframe for all taxa with occurrences >= 10:
  streptophyta_sib55_occ10 <- streptophyta_sib55[streptophyta_sib55$scientificName %in% gbif_sib55_occurrences10$scientificName, ]
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### Taxa occurrences for "Beringia"-Region (40-90°N, 150-220°E):

{
  
  # read occurrence table:
  gbif_beringia <- read.table("input/gbif_taxa_occurences_Beringia_region/occurrence.txt", 
                              sep="\t", header=TRUE, comment.char="#", na.strings=".", stringsAsFactors=FALSE, quote="", fill=TRUE)
  
  
  # filter the gbif_beringia dataframe for only those columns that should be kept:
  gbif_beringia_keep <- gbif_beringia[ ,which(names(gbif_beringia) %in% keep)]
  
  
  # filter taxa in the gbif_beringia_keep dataframe for Streptophyta:
  streptophyta_beringia <- subset(gbif_beringia_keep, phylum == "Tracheophyta" | phylum == "Bryophyta" | 
                                    phylum == "Marchantiophyta" | phylum == "Anthocerotophyta")
  
  # -------------------------------------------------------------------------------------------------
  
  # calculate number of occurrences for each taxon: 
  
  beringia_IDs <- sort(unique(streptophyta_beringia$scientificName))
  
  gbif_beringia_occurrences <- data.frame(acceptedScientificName=NA, scientificName=NA, family=NA, genus=NA, species=NA, occurrence=NA)
  
  
  for(i in 1:length(beringia_IDs)){
    
    print(paste0(i,"/",length(beringia_IDs)))
    
    taxasubset <- streptophyta_beringia[streptophyta_beringia$scientificName == beringia_IDs[i], ]
    
    accScienNam <- unique(taxasubset$acceptedScientificName)
    spec        <- unique(taxasubset$species)
    
    gbif_beringia_occurrences[i,"acceptedScientificName"] <- accScienNam[which(accScienNam != "")]
    gbif_beringia_occurrences[i,"scientificName"]         <- unique(taxasubset$scientificName)
    gbif_beringia_occurrences[i,"family"]                 <- unique(taxasubset$family) 
    gbif_beringia_occurrences[i,"genus"]                  <- unique(taxasubset$genus)
    gbif_beringia_occurrences[i,"species"]                <- ifelse(length(spec) > 1, spec[which(spec != "")], spec)
    gbif_beringia_occurrences[i,"occurrence"]             <- dim(taxasubset)[1]
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # filter occurrence dataframe for all taxa with occurrences >= 10:
  gbif_beringia_occurrences10 <- subset(gbif_beringia_occurrences, occurrence >= 10)
  
  # filter streptophyta_beringia dataframe for all taxa with occurrences >= 10:
  streptophyta_beringia_occ10 <- streptophyta_beringia[streptophyta_beringia$scientificName %in% gbif_beringia_occurrences10$scientificName, ]
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### combine "Siberia55" and "Beringia" dataframes:

{
  
  # combine "Siberia55" and "Beringia" dataframes:
  streptophyta_combined <- rbind(streptophyta_sib55, streptophyta_beringia)
  
  
  # calculate number of occurrences for each taxon: 
  
  combined_IDs <- sort(unique(streptophyta_combined$scientificName))
  
  gbif_combined_occurrences <- data.frame(acceptedScientificName=NA, scientificName=NA, family=NA, genus=NA, species=NA, occurrence=NA)
  
  
  for(i in 1:length(combined_IDs)){
    
    print(paste0(i,"/",length(combined_IDs)))
    
    taxasubset <- streptophyta_combined[streptophyta_combined$scientificName == combined_IDs[i], ]
    
    accScienNam <- unique(taxasubset$acceptedScientificName)
    spec        <- unique(taxasubset$species)
    
    gbif_combined_occurrences[i,"acceptedScientificName"] <- accScienNam[which(accScienNam != "")]
    gbif_combined_occurrences[i,"scientificName"]         <- unique(taxasubset$scientificName)
    gbif_combined_occurrences[i,"family"]                 <- unique(taxasubset$family) 
    gbif_combined_occurrences[i,"genus"]                  <- unique(taxasubset$genus)
    gbif_combined_occurrences[i,"species"]                <- ifelse(length(spec) > 1, spec[which(spec != "")], spec)
    gbif_combined_occurrences[i,"occurrence"]             <- dim(taxasubset)[1]
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # filter occurrence dataframe for all taxa with occurrences >= 10:
  gbif_combined_occurrences10 <- subset(gbif_combined_occurrences, occurrence >= 10)
  
  # filter streptophyta_combined dataframe for all taxa with occurrences >= 10:
  streptophyta_combined_occ10 <- streptophyta_combined[streptophyta_combined$scientificName %in% gbif_combined_occurrences10$scientificName, ]
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

save(streptophyta_combined, file = "output/steptophyta_combined.RData")
