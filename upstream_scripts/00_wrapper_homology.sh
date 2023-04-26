#!/bin/bash

# MolEvolvR: Companion wrapper script
# accnum –> fasta –> deltablast -> edirect -> blastclust -> iprscan/rpsblast -> cleanup w/ lineages!
# To run on compute.cvm.msu.edu

# Created: 2020.07.09
# Last modified: 2020.12.14
# Authors: Lauren Sosinski, Janani Ravi
# On GitHub: currently in jravilab/molevol_scripts/upstream_scripts

# Input type: text file containing paths to individual fasta files to analyze

############
## TORQUE ##
############
#PBS -l nodes=1:ppn=10		# number of nodes requested, specify wall time
#PBS -m abe			# email notifications for job
#PBS -M=jravilab.msu@gmail.com	# user email;
#PBS -N molevol_analysis	# name of job being run

## print start/stop printf in individual scripts

## change output directory based on user input
OUTPATH=$PBS_O_WORKDIR
cd ${OUTPATH}

## !! NOTES FOR LS ##
## add flags/help function after finishing script and/or gencontext
## talk to Sam about how this works with Pins package

## USER INPUTS
INPUTPATHS_LIST=$1
DB=$2
NHITS=$3
EVAL=$4
IS_QUERY=$5


## USAGE
## qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F "input.txt F"
## qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F "example_blastp.csv T"

# Location of databases/dependencies
export BLASTDB=/data/common_data/blastdb/v6
export BLASTMAT=/opt/software/BLAST/2.2.26/data
export INTERPRO=/opt/software/iprscan/5.47.82.0-Python3/data:/data/common_data/iprscan:/var/interproscan/data:$INTERPRO
export SIGNALP=/var/interproscan/bin/signalp/4.1
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

#total=$(( ((3*60)+50)*60 ))			## total time the job can take in seconds, this should match your SBATCH line above
#maxtime=$(( 120*60 )) 				## maximum time to process one input, need to do some experimenting with your inputs
if [ "$IS_QUERY" = "T" ]
# Handle the query data
then
FILE="$INPUTPATHS_LIST"
F=$(basename ${FILE})
PREFIX="query_data"
OUTDIR=${OUTPATH}/${PREFIX}
mkdir ${OUTDIR}				## make the directory
cat $FILE >> ${OUTDIR}/${PREFIX}.all_accnums.fa
if [ -e starting_accs.txt ]
then
cp starting_accs.txt ${OUTDIR}/${PREFIX}.all_accnums.txt
else
cp accs.txt ${OUTDIR}/${PREFIX}.all_accnums.txt
fi
cd ${OUTDIR} || exit
#Acc2info
sh /data/research/jravilab/molevol_scripts/upstream_scripts/acc2info.sh ${OUTDIR}/${PREFIX}.all_accnums.txt $PREFIX $OUTDIR
mv ${OUTDIR}/${PREFIX}.acc2info.tsv ${OUTDIR}/${PREFIX}.blast.cln.tsv
else
# Handle homolog data
	FILE=$(sed -n "${PBS_ARRAYID}"p "${INPUTPATHS_LIST}")
	F=$(basename ${FILE})
   	PREFIX=$(echo "${F%%.faa}")	## takes PREFIX of file
   	OUTDIR=${OUTPATH}/${PREFIX}_full	## variable containing output filepath based PREFIX
   	printf "${PREFIX}\n"
	mkdir ${OUTDIR}				## make the directory
    cd ${OUTDIR} || exit

	## DELTABLAST ##
	db_start=$SECONDS
	sh /data/research/jravilab/molevol_scripts/upstream_scripts/01.1_deltablast.sh $FILE $PREFIX $OUTDIR $DB $NHITS $EVAL
	db_dur=$(( $SECONDS - $db_start ))

    ## ACC2FA -- getting fasta FILES for deltablast output(s)
	acc2fa_start=$SECONDS
	sh /data/research/jravilab/molevol_scripts/upstream_scripts/02_acc2fa.sh ${OUTDIR}/${PREFIX}.dblast.tsv $PREFIX $OUTDIR
	acc2fa_dur=$(( $SECONDS - $acc2fa_start ))

	## ACC2INFO ##
	acc2info_start=$SECONDS
	sh /data/research/jravilab/molevol_scripts/upstream_scripts/acc2info.sh ${OUTDIR}/${PREFIX}.all_accnums.txt $PREFIX $OUTDIR
	acc2info_dur=$(( $SECONDS - $acc2info_start ))

	## BLAST RESULT CLEANUP ##
	db_cln_start=$SECONDS
	Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/01.2_cleanup_blast.R ${OUTDIR}/${PREFIX}.dblast.tsv ${OUTDIR}/${PREFIX}.acc2info.tsv $PREFIX $F
	db_cln_dur=$(( $SECONDS - $db_cln_start ))
fi 
	## BLASTCLUST ##
	bclust_start=$SECONDS
	sh /data/research/jravilab/molevol_scripts/upstream_scripts/03.1_blastclust.sh ${OUTDIR}/${PREFIX}.all_accnums.fa $PREFIX $OUTDIR
	bclust_dur=$(( $SECONDS - $bclust_start ))

	## CLUST2TABLE
	c2t_start=$SECONDS
	Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/03.2_clust2table.R ${OUTDIR}/${PREFIX}.bclust.L60S80.tsv ${OUTDIR}/${PREFIX}.blast.cln.tsv
	c2t_dur=$(( $SECONDS - $c2t_start ))

	## INTERPROSCAN ##
	## add second run for original protein, too
	ipr_start=$SECONDS
	sh /data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh ${OUTDIR}/${PREFIX}.all_accnums.fa ${PREFIX} ${OUTDIR}
	ipr_dur=$(( $SECONDS - $ipr_start ))

	## IPR2LIN ##
	#Append colnames to beginning of ipr file
	sed -i '1s/^/AccNum\tSeqMD5Digest\tSLength\tAnalysis\tDB.ID\tSignDesc\tStartLoc\tStopLoc\tScore\tStatus\tRunDate\tIPRAcc\tIPRDesc\n/' ${OUTDIR}/${PREFIX}.iprscan.tsv
	ipr2lin_start=$SECONDS
	Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/01.4_ipr2lin.R ${OUTDIR}/${PREFIX}.iprscan.tsv ${OUTDIR}/${PREFIX}.acc2info.tsv $PREFIX
	ipr2lin_dur=$(( $SECONDS - $ipr2lin_start ))

	## IPR2DA ##
	ipr2da_start=$SECONDS
	Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/05a_ipr2da.R ${OUTDIR}/${PREFIX}.iprscan_cln.tsv ${PREFIX} ${OUTDIR}/${PREFIX}.cln.clust.tsv
	ipr2da_dur=$(( $SECONDS - $ipr2da_start ))

	## RPSBLAST ##
	rps_start=$SECONDS
	#sh /data/research/jravilab/molevol_scripts/upstream_scripts/04b_rpsblast.sh ${OUTDIR}/${PREFIX}.all_accnums.fa ${PREFIX} ${OUTDIR}
	rps_dur=$(( $SECONDS - $rps_start))

	## RPS2DA ##
	rps2da_start=$SECONDS
	#Rscript /data/research/jravilab/molevol_scripts/upstream_scripts/05b_rps2da.R ${OUTDIR}/${PREFIX}.rps.out ${OUTDIR}/${PREFIX}.cln.clust.ipr.tsv ${PREFIX}
	rps2da_dur=$(( $SECONDS - $rps2da_start ))

	cp ${FILE} ${OUTDIR}	# copy fasta file to output directory

	## Figure out how long the entire script took to run
	dur=$(( $SECONDS - $start ))
	printf "\nTotal run time: $dur\n"
	STOP_DT=$(date '+%d/%m/%Y-%H:%M:%S')

	## Add benchmarking times to logfile
	printf "\n${START_DT}\t${STOP_DT}\t${PREFIX}\t${db_dur}\t${acc2info_dur}\t${db_cln_dur}\t${acc2fa_dur}\t${bclust_dur}\t${c2t_dur}\t${ipr_dur}\t${ipr2da_dur}\t${dur}" >> ${OUTPATH}/logfile.tsv

	## And how much time is left
	#timeleft=$(( $total - $dur ))

        ## If there's a chance we get a long input to process, then
        ## resubmit this job, then kill this job
        #     if [ ${timeleft} -lt ${maxtime} ]; then
        #       sbatch --array=${SLURM_ARRAY_TASK_ID} $0
        #      scancel ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
        #  fi
NUM_RUNS=$(wc -l "${OUTPATH}"/logfile.tsv | grep -Eo "^[[:digit:]]+")
((NUM_RUNS-=1))
TOTAL_RUNS=$(wc -l "${OUTPATH}"/input.txt | grep -Eo "^[[:digit:]]+")
if [ $TOTAL_RUNS -eq $NUM_RUNS ]
then
  touch ../done.txt
else
  echo "${NUM_RUNS} / ${TOTAL_RUNS} jobs completed" > ../status.txt
fi


# FA: disabled b/c docker's handling perms
# setfacl -R -m group:shiny:r-x ${OUTDIR}