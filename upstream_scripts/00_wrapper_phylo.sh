#!/bin/bash
# MolEvolvR: Companion wrapper script
# accnum –> fasta –> iprscan/rpsblast -> cleanup w/ lineages!
# To run on compute.cvm.msu.edu
# Created: 2021.12.12
# Last modified: 2021.12.12
# Authors: Joe Burke
# On GitHub: currently in jravilab/molevol_scripts/upstream_scripts
# Input type: text file containing paths to individual fasta files to analyze
############
## TORQUE ##
############
#PBS -l nodes=1:ppn=10,walltime=10:00:00		# number of nodes requested, specify wall time
#PBS -m abe			# email notifications for job
#PBS -M=jravilab.msu@gmail.com	# user email;
#PBS -N molevol_analysis_phylo	# name of job being run
# Begin
OUTPATH=$PBS_O_WORKDIR
cd ${OUTPATH}

INPUTPATHS_LIST=$1
DB=refseq
NHITS=1000
EVAL=1e-5
WBLAST=$2

# Location of databases/dependencies
export BLASTDB=/data/common_data/blastdb/v5:/data/common_data/blastdb/ncbidb:/data/common_data/blastdb:/data/common_data/blastdb/FASTA:$BLASTDB
export BLASTMAT=/opt/software/BLAST/2.2.26/data
export INTERPRO=/opt/software/iprscan/5.47.82.0-Python3/data:/data/common_data/iprscan:$INTERPRO
export NCBI_API_KEY=882b28aa19ece4679d4fa5adcf3319f5df09


#####################
## LOADING MODULES ##
#####################

# Prevent "module: command not found"
# Read more about it https://www.sdsc.edu/support/user_guides/tscc.html
source /etc/profile.d/modules.sh

module purge 					## clear loaded modules
module load R		 			## load R
module load edirect 				## load edirect
module load BLAST				## load blast (for blastclust)
module load BLAST+ BioPerl 			## load blast+
module load iprscan 				## load iprscan

start=$SECONDS 					## get current time
START_DT=$(date '+%d/%m/%Y-%H:%M:%S')

FILE=$(sed -n "${PBS_ARRAYID}"p "${INPUTPATHS_LIST}")

F=$(basename "${FILE}")
PREFIX="${F%%.faa}"		## takes PREFIX of file
OUTDIR=${OUTPATH}/${PREFIX}_phylo
mkdir ${OUTDIR}		## variable containing output filepath based PREFIX
cd "${OUTDIR}" || exit

cp ${FILE} ${OUTDIR}	# copy fasta file to output directory
## ACC2INFO ##
acc2info_start=$SECONDS
printf "AccNum.noV\tFullAccNum\tDescription\tLength\tTaxID\tSpecies\tSourceDB\tCompleteness\n" > ${OUTDIR}/${PREFIX}.acc2info.tsv
grep "${PREFIX}" ../acc2info.tsv >> ${OUTDIR}/${PREFIX}.acc2info.tsv
acc2info_dur=$(( $SECONDS - $acc2info_start ))

## INTERPROSCAN ##
## add second run for original protein, too
ipr_start=$SECONDS
sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh ${OUTDIR}/${PREFIX}.faa ${PREFIX} ${OUTDIR}
ipr_dur=$(( $SECONDS - $ipr_start ))

## IPR2LIN ##
ipr2lin_start=$SECONDS
Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/01.4_ipr2lin.R ${OUTDIR}/${PREFIX}.iprscan.tsv ${OUTDIR}/${PREFIX}.acc2info.tsv $PREFIX
ipr2lin_dur=$(( $SECONDS - $ipr2lin_start ))

## IPR2DA ##
ipr2da_start=$SECONDS
Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/05a_ipr2da.R ${OUTDIR}/${PREFIX}.iprscan_cln.tsv ${PREFIX} ${OUTDIR}/${PREFIX}.acc2info.tsv
ipr2da_dur=$(( $SECONDS - $ipr2da_start ))

## Figure out how long the entire script took to run
dur=$(( $SECONDS - $start ))
printf "\nTotal run time: $dur\n"
STOP_DT=$(date '+%d/%m/%Y-%H:%M:%S')


## Add benchmarking times to logfile
printf "\n${START_DT}\t${STOP_DT}\t${PREFIX}\t${acc2info_dur}}\t${ipr_dur}\t${ipr2da_dur}\t${dur}" >> ${OUTPATH}/logfile.tsv

NUM_RUNS=$(wc -l "${OUTPATH}"/logfile.tsv | grep -Eo "^[[:digit:]]+")
((NUM_RUNS-=1))
TOTAL_RUNS=$(wc -l "${OUTPATH}"/input.txt | grep -Eo "^[[:digit:]]+")
if [ "$TOTAL_RUNS" -eq $NUM_RUNS ]
then
  touch ../done.txt
else
  echo "${NUM_RUNS} / ${TOTAL_RUNS} jobs completed" > ../status.txt
fi