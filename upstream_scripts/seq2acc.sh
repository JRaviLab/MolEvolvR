QUERY_FA=$1
# Location of databases/dependencies 
export BLASTDB=/data/common_data/blastdb/v6
export BLASTMAT=/opt/software/BLAST/2.2.26/data

source /etc/profile.d/modules.sh 
module purge 				
module load BLAST+				

OUTPATH="$PBS_O_WORKDIR""/seq2acc_data"
INPUT_FILE_NAME=$(basename "$QUERY_FA" .fa)

blastp -query "$QUERY_FA" -db refseq_protein -out "$OUTPATH"/"$INPUT_FILE_NAME".csv \
    -evalue 0.000001 -max_target_seqs 5 -num_threads 20 
#$PC_IDENTITY=$(awk -F "," 'NR<2 { print $3 }' test.csv)
#pident is 3rd column https://www.ncbi.nlm.nih.gov/books/NBK279684/
