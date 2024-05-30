# sedaDNA-g-h-metabarcoding-SibAla_2023

## SibAla_2023 database creation
Scripts 0 to 2 need to be executed one after the other.

0_GBIF_taxa_from_the_area.R: used to compile all species of Streptophyta with a minimum of 10 occurences in Alaska and Siberia (55-90°N, 50-150°E ; 40-90°N, 150°E-140°W).
Need the downloaded occurence tables:
- 2022-07-08_sequencetable_taxanames_lineage_siberia55_occ10_final.csv
- gbif_taxa_occurences_sib55_region.txt
- gbif_taxa_occurences_Beringia_region.txt

1_GBIF_to_NCBI_taxonomy.R: will make compatible both GBIF and NCBI taxonomis to keep only species occuring in our region and in the reference databases (arctborbryo, EMBL143 and PhyloNorway).
Need:
- output from previous script
- taxdump taxonomy downloaded from NCBI
- arctoborbryo database prepared for OBITools (trnL g/h): 2022-05-05_arcbryo_GH_datatable.csv
- phylonorway database prepared for OBITools (trnL g/h): 2022-05-05_phylonorway_GH_datatable.csv
- embl143 database prepared for OBITools (trnL g/h): 2022-05-11_embl143_GH_datatable.csv
- List of contaminants species to remove: KS_conta_list.csv

2_data_to_fasta.R: clean the previous datatable and convert to fasta to provide as input to create an OBITools3 database.
  - Only the out from the previous script is needed.

## OBITools3 against SibAla_2023 


## Quality filtering of the replicates


