# sedaDNA-g-h-metabarcoding-SibAla_2023

## SibAla_2023 database creation
Scripts 0 to 2 need to be executed one after the other.

0_GBIF_taxa_from_the_area.R: used to compile all species of Streptophyta with a minimum of 10 occurrences in Alaska and Siberia (55-90°N, 50-150°E ; 40-90°N, 150°E-140°W) from the  Global Biodiversity Information Facility database (GBIF, ).
Need the downloaded occurence tables:
- 2022-07-08_sequencetable_taxanames_lineage_siberia55_occ10_final.csv
- gbif_taxa_occurences_sib55_region.txt
- gbif_taxa_occurences_Beringia_region.txt

1_GBIF_to_NCBI_taxonomy.R: will make compatible both GBIF and NCBI taxonomies to keep only species occurring in the study region of Siberia and Alaska and in the reference databases (arctborbryo (2,3,4), EMBL143 (5) and PhyloNorway (6)).
Need:
- output from previous script
- taxdump taxonomy downloaded from NCBI
- arctborbryo database prepared for OBITools (trnL g/h): 2022-05-05_arcbryo_GH_datatable.csv
- PhyloNorway database prepared for OBITools (trnL g/h): 2022-05-05_phylonorway_GH_datatable.csv
- embl143 database prepared for OBITools (trnL g/h): 2022-05-11_embl143_GH_datatable.csv
- List of contaminants species to remove: KS_conta_list.csv

2_data_to_fasta.R: clean the previous data table and convert to fasta to provide as input to create an OBITools3 database.
- Only the output from the previous script is needed.

3_SibAla_2023_build.sl: build the SibAla_2023 database for OBITools.
- using the fasta file output from previous script.

## OBITools3 against SibAla_2023 


## Quality filtering of the replicates

## References
1. GBIF.org (2023), GBIF Home Page. Available from: https://www.gbif.org [30 May 2024].
2. E. Willerslev, J. Davison, M. Moora, M. Zobel, E. Coissac, M. E. Edwards, E. D. Lorenzen, M. Vestergård, G. Gussarova, J. Haile, Fifty thousand years of Arctic vegetation and megafaunal diet. Nature, 506, 47-51 (2014). https://doi.org/10.1038/nature12921
3. E. M. Soininen, G. Gauthier, F. Bilodeau, D. Berteaux, L. Gielly, P. Taberlet, G. Gussarova, E. Bellemain, K. Hassel, H. K. Stenøien, Highly overlapping winter diet in two sympatric lemming species revealed by DNA metabarcoding. PLoS One, 10, e0115335 (2015). https://doi.org/10.1371/journal.pone.0115335 
4. J. Sønstebø, L. Gielly, A. Brysting, R. Elven, M. Edwards, J. Haile, E. Willerslev, E. Coissac, D. Rioux, J. Sannier, Using next‐generation sequencing for molecular reconstruction of past Arctic vegetation and climate. Molecular Ecology Resources, 10, 1009-1018 (2010). https://doi.org/10.1111/j.1755-0998.2010.02855.x
5. C. Kanz, P. Aldebert, N. Althorpe, W. Baker, A. Baldwin, K. Bates, P. Browne, A. van den Broek, M. Castro, G. Cochrane, The EMBL nucleotide sequence database. Nucleic Acids Research 33, D29-D33 (2005). https://doi.org/10.1093/nar/gki098
6. I. G. Alsos, D. P. Rijal, D. Ehrich, D. N. Karger, N. G. Yoccoz, P. D. Heintzman, A. G. Brown, Y. Lammers, L. Pellissier, T. Alm, Postglacial species arrival and diversity buildup of northern ecosystems took millennia. Science Advances, 8, eabo7434 (2022). https://doi.org/10.1126/sciadv.abo7434
