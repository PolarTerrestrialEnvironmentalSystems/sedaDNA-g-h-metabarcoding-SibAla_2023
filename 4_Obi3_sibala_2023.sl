#install and load obitools3
module load obitools/3.0.1b20

#data variables *enter name of the sequencing run ALRK-5, ALRK-11 etc.
DMS_OUT=output/out
DMS=APMG-68
DATAPATH=datapath/data
FORWARD=*_R1.fastq
REVERSE=*_R2.fastq
TAGFILE=$DATAPATH/*_tagfile.txt
ID=*
#reference database sibala_2023
SIBALA=/albedo/work/projects/p_biodiv_dbs/obitools/obi.obidms

# set tmp variables
DATA1=$DATAPATH/$FORWARD
DATA2=$DATAPATH/$REVERSE

## import data
srun obi import --fastq-input $DATA1 $DMS/${ID}reads1
srun obi import --fastq-input $DATA2 $DMS/${ID}reads2
## import tag file
srun obi import --ngsfilter-input $TAGFILE $DMS/${ID}tagfile
## align paired ends
srun obi alignpairedend -R $DMS/${ID}reads1 $DMS/${ID}reads2 $DMS/${ID}aligned_reads
## remove unaligned
srun obi grep -a mode:alignment $DMS/${ID}aligned_reads $DMS/${ID}good_sequences
## assign reads to tag combinations
srun obi ngsfilter -t $DMS/${ID}tagfile -u $DMS/${ID}unidentified_seqs $DMS/${ID}good_sequences $DMS/${ID}identified_sequences
# dereplicate sequences into unique sequences
srun obi uniq --merge sample $DMS/${ID}identified_sequences $DMS/${ID}dereplicated_sequences
## denoise data, only keep COUNT and merged_sample tags
srun obi annotate -k COUNT -k MERGED_sample $DMS/${ID}dereplicated_sequences $DMS/${ID}cleaned_metadata_sequences
## grep
srun obi grep -p "len(sequence)>=1 and sequence['COUNT']>=1" $DMS/${ID}cleaned_metadata_sequences $DMS/${ID}cleaned_10_metadata_sequences
## denoise data, clean from pcr/sequencing errors
srun obi clean -s MERGED_sample -r 0.05 -H $DMS/${ID}cleaned_10_metadata_sequences $DMS/${ID}cleaned_sequences

### SIBALA_2023 reference database

# clean DMS
srun obi clean_dms $DMS
# assign to SIB_db
srun obi clean_dms $SIBALA
srun obi ecotag --taxonomy $SIBALA/TAXONOMY/ncbi_tax -R $SIBALA/VIEWS/2023-03-15_sibala_97_DB.obiview $DMS/${ID}cleaned_sequences $DMS/${ID}sib23_assigned_sequences
# annonate lineage information
srun obi annotate --with-taxon-at-rank family --with-taxon-at-rank genus --with-taxon-at-rank species --taxonomy $SIBALA/TAXONOMY/ncbi_tax $DMS/${ID}sib23_assigned_sequences $DMS/${ID}sib23_annotated_sequences
# export embl assignment as csv file
srun obi export --tab-output $DMS/${ID}sib23_annotated_sequences > ${DMS}.obidms/${ID}sib23_anno.csv

# copy dms (all files of analysis) from calculating node to work folder
srun mv ${DMS}.obidms $DMS_OUT

