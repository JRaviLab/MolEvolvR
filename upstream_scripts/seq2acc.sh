QUERY_FA=$1
# Location of databases/dependencies 
export BLASTDB=/data/common_data/blastdb/v6
export BLASTMAT=/opt/software/BLAST/2.2.26/data
export INTERPRO=/opt/software/iprscan/5.47.82.0-Python3/data:/data/common_data/iprscan:$INTERPRO
export NCBI_API_KEY=YOUR_KEY_HERE

#####################
## LOADING MODULES ##
#####################

# Prevent "module: command not found"
# Read more about it https://www.sdsc.edu/support/user_guides/tscc.html
source /etc/profile.d/modules.sh 
module purge 					## clear loaded modules
module load BLAST+				## load blast (for blastclust				

OUTPATH="$PBS_O_WORKDIR""/seq2acc_data"
INPUT_FILE_NAME=$(basename "$QUERY_FA" .fa)

blastp -query "$QUERY_FA" -db refseq_protein -out "$OUTPATH"/"$INPUT_FILE_NAME".csv \
    -evalue 0.000001 -max_target_seqs 5 -num_threads 6 
#$PC_IDENTITY=$(awk -F "," 'NR<2 { print $3 }' test.csv)
#pident is 3rd column https://www.ncbi.nlm.nih.gov/books/NBK279684/
