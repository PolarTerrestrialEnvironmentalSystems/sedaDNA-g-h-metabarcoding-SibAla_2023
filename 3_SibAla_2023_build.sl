# load obitools3
module load bio/OBItools/3.0.1b16

# build obitools sibala_2023 database
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
obi import --taxdump taxdump.tar.gz obi/taxonomy/ncbi_tax
obi import --fasta 2023-03-17_Sib5510Seq_sequencetable_siberia55_occ10_clean_nodots_withrank_final.fasta obi/2023-03-15_Sibala
obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy obi/taxonomy/ncbi_tax obi/2023-03-15_sibala obi/2023-03-15_sibala_clean
obi uniq --taxonomy obi/taxonomy/ncbi_tax obi/2023-03-15_sibala_clean obi/2023-03-15_sibala_uniq
obi build_ref_db -t 0.97 --taxonomy obi/taxonomy/ncbi_tax obi/2023-03-15_sibala_uniq obi/2023-03-15_sibala_97_DB
