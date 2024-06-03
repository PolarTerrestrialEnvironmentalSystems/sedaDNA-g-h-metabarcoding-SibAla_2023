
### R-Script for filtering and resampling (on ASV, best taxonomic name or family level) the sibala_2023 data files 

rm(list=ls())
library(readxl)
library(stringr)
library(tidyverse)

# helping function:
`%notin%` <- Negate(`%in%`)

### HEADER preparation -  define your in- and output paths 
  
inputpath  <- "define your path"
outputpath90 <- "define your path"
# read in the according agefile for the dataset 
age_depth  <- read_xlsx("*_agefile.xlsx", sheet="Sheet1", col_names=TRUE, trim_ws=TRUE)
# read in a list of PCR replicates with bad quality to remove the samples from the dataset before resampling
replicates_remove <- read_xlsx("*_List_rep_cont_to_remove.xlsx", sheet="Sheet1", col_names=TRUE, trim_ws=TRUE)


  # select sequencing run (ALRK-5, ALRK-10, HUA-9 etc.):
  corename <-  "HUA-9"
  
  # select lake core:
  lake <- "PG2133"
  
  # select identity level for which you want to filter the data:
  identity_level <- 0.9
  

  # remove families that contain "NA": 
  ### please select "yes" or "no"
  remove_NAfamilies <- "yes"
  

  # save output files: 
  ### please select "yes" or "no"
  save_output <- "yes"
  

### PART 1: read in data files from Obitools 3 pipeline
  

    # read data:
    arc      <- read.csv(paste0(inputpath,corename,"_sibala2023_anno.csv"), header=TRUE, sep="\t", dec=".", check.names="FALSE")
    arc_only <- arc %>% select(starts_with("MERGED_sample:"),NUC_SEQ, ID, BEST_IDENTITY, TAXID, SCIENTIFIC_NAME, family, family_name, genus, genus_name, species, species_name)
    
    arc_only_df <-arc_only %>% 
      rename(
        NUC_SEQ_arc=NUC_SEQ,
        best_identity=BEST_IDENTITY,
        scientific_name=SCIENTIFIC_NAME,
        best_family=family_name,
        best_genus=genus_name,
        best_species=species_name
      )
    names(arc_only_df)
    arc_only_df$best_identity
    
    dim(arc_only_df)

    arc_reprem_df <- arc_only_df[ ,names(arc_only_df) %notin% replicates_remove$sample_remove]      
    str(arc_reprem_df$best_identity)
    dim(arc_reprem_df)

    # filter dataset for best identity = 90%:
    arc_embl_combi_filter <- arc_reprem_df %>% filter(best_identity >= identity_level)
   
    # select only columns needed:
    arc_embl_combi_filter <- arc_embl_combi_filter %>% select(best_identity, NUC_SEQ_arc, scientific_name, best_family, best_genus, best_species,
                                                              starts_with("MERGED_sample:"))
    
    # filter dataset for lake:
    arc_embl_combi_filter_lake <- cbind(arc_embl_combi_filter[ ,c(1:6)], arc_embl_combi_filter[, (grepl(lake, colnames(arc_embl_combi_filter)))])


### PART 2: preparation for resampling:


  # remove families that are "NA":
  if(remove_NAfamilies == "yes"){
    
    arc_embl_combi_filter_lake <- arc_embl_combi_filter_lake[!is.na(arc_embl_combi_filter_lake$best_family), ]
    
  }

  names(arc_embl_combi_filter_lake)
  sampleid <- names(select(arc_embl_combi_filter_lake, contains("sample:")))
  
  # extract sequence, family, samples and blanks:
  cols_need <- c("NUC_SEQ_arc", "scientific_name", "best_family", sampleid) 
  id1_need  <- arc_embl_combi_filter_lake[, cols_need] 
  
  
  # select only blanks (for having a look on only the blanks):
  id1_blanks <- data.frame(id1_need[ ,c(1:3)], id1_need[, (grepl("blank", colnames(id1_need)))])
  id1_blanks[ ,-c(1:3)][is.na(id1_blanks[ ,-c(1:3)])] <- 0
  
  # remove blanks:
  id1_samples <- id1_need[, !(grepl("blank", colnames(id1_need)))]
  
  # set NA == 0:
  id1_count <- id1_samples[-c(1:3)]
  id1_count[is.na(id1_count)] <- 0
  
  # extract the Extraction_number and set them as column names:
  colna_extract <- unlist(lapply(strsplit(split="_", names(id1_count)),function(x)return(x[4])))  # the extraction number is the forth element after splitting up the sample name
  colnames(id1_count) <- colna_extract
  
  # make the column names unique (for identifying replicates):
  colnames(id1_count) <- make.unique(as.character(names(id1_count)), sep = "_")
  names(id1_count)
  
  # sort column names:
  id1_count_sort <- id1_count %>% select(sort(names(.)))
  
  # attach sequence and scientific_name:
  id1_sort <- cbind(id1_samples[c(1:3)], id1_count_sort)
  
  # remove sequences with total count = 0:
  id1_sort_done <- id1_sort[apply(id1_sort[ ,-c(1:3)], 1, function(x) !all(x==0)), ]
  
  # remove sequences with scientific_name == NA or best_family == NA:
  final <- id1_sort_done[complete.cases(id1_sort_done[ ,c(2,3)]), ] 
  
  # check if there are still NA in the dataframe (needs to be "integer(0)" if everything is ok):
  which(is.na(final))
  

  # merge the replicates by extraction number:
  
  # find the replicates:
  replicates <- unlist(lapply(strsplit(split="_", names(final)[-c(1:3)]),function(x)return(x[1])))
  names(final)[-c(1:3)] <- replicates
  names(final)
  
  # extract the unique extraction number:
  samplespresent <- unique(replicates)
  
  # merge:
  specseq <- final
  specseq_merge <- final[c(1:3)]
  
  for(samplei in samplespresent){
    
    print(samplei)
    
    clnumberwithrepeat <- grep(x=names(specseq), pattern=samplei)
    
    if(length(clnumberwithrepeat) > 0){
      
      if(length(clnumberwithrepeat) > 1){
        
        dfi <- data.frame(apply(specseq[ ,clnumberwithrepeat], 1, sum))
        names(dfi) <- samplei
        specseq_merge <- cbind(specseq_merge, dfi)
        
      }else{
        
        dfi <- data.frame(specseq[ ,which(names(specseq) == samplei)])
        names(dfi) <- samplei
        specseq_merge <- cbind(specseq_merge, dfi)
      }
      
    }
    
  }
  
  str(specseq_merge)
  names(specseq_merge)
  
  # save output table:
  if(save_output == "yes"){write.csv2(specseq_merge, file=paste0(outputpath90,corename,"_",lake,"_file02_identitylevel_",identity_level*100,"_merged_replicates.csv"), row.names=FALSE)}

  ### add better name with depth and age and add taxonomic information from family and scientific name
  
  # connect data frame with age-depth information, here starts_with("ESB") needs to be adjusted in each dataset
  specseq_merge_long=pivot_longer(data=specseq_merge,
                                  cols = starts_with("ESB"),
                                  names_to = "Extraction_number",
                                  values_to = "counts")
  
  specseq_merge_long_ad <- left_join(specseq_merge_long, unique(age_depth[ ,c(3,6,7)]))
  
  # merge extraction number, depth and age into a new sample_name and add to data frame:
  specseq_merge_long_ad_colmerge <- unite(specseq_merge_long_ad, Extraction_number, Mean_collection_depth_round, Section_age_round, col = "sample_name", sep = "_")
  names(specseq_merge_long_ad_colmerge)
  
  specseq_merge_long_age_wide <- pivot_wider(data=specseq_merge_long_ad_colmerge[ ,c(1,4,5)], 
                                             names_from = "sample_name",  
                                             values_from = "counts")
  
  specseq_merge_long_age_melt=full_join(specseq_merge_long_age_wide, arc_embl_combi_filter_lake[ ,c(1,2,3,4,5,6)], by="NUC_SEQ_arc")
  
  # save output table:
  # sample names are extraction numbers and depth and ages can be added via age-depth file from the "metadata folder":
  if(save_output == "yes"){write.csv2(specseq_merge_long_age_melt, file=paste0(outputpath90,corename,"_",lake,"_file02.2_identitylevel_",identity_level*100,"_merged_replicates.csv"), row.names=FALSE)}
