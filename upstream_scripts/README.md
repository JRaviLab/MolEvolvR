# Scripts and output files

## Scripts

### Homology search
#### DELTABLAST
`dblast.sb`/`deltablast.sb`/`deltablast.sh` <br>
- input: fasta file(s)
- output: file containing similar proteins based on sequence similarity; `org_prot.DB.1e5.out`
- columns: qseqid sgi sseqid pident length mismatch qstart qend sstart send evalue bitscore staxids
- database (DB): nr = non-redundant sequences, refseq = reference sequences
- evalue: 1e^(-5)
- does not contain arguments, call script with `sbatch dblast.sb INFILE`

#### RPSBLAST
`rpsblast.sb`/`rpsblast.sh` <br>
- performs RPSBlast analysis on file containing fasta sequences from CDD database
- input: file containing fasta sequences; `org_prot.all_accnums.fa`
- output: similar to deltablast output, but contains CDD domains for each protein; `org_prot.EVAL.rps.txt`

#### rps2viz

#### rps2da
`rps2da.R` <br>
- input: rpsblast output, deltablast output
- output: blast output w/ domain architecture added as a string
- still needs a couple of bugs fixed + eggNOG database implemented (instead of incomplete cdd)

### acc2fa
`acc2fa.sh` <br>
- input: txt file containing deltablast results for protein `org_prot.DB.EVAL.txt`
- output: fasta sequences in single file to be pushed through clustering program
- uses entrez direct to grab fasta files for each homolog

### Clustering
#### BLASTCLUST
`blastclust.sb`/`blastclust.sh` <br>
- input: file of fasta sequences, given on command line `org_prot.all_accnums.fa`
- output: file with clusters of similar proteins `org_prot.bclust.LnSm.out`
- produces outputs for L60S80 and L90S95 clustering parameters

#### BLAST+BLASTCLUST
`clust2table.R` <br>
- script to generate cluster names based on position in `blastclust` output and how many proteins are present in the cluster
- input: `blastclust` results, `deltablast` results; `org_prot.DB.EVAL.txt` and `org_prot.bclust.LnSm.out`
- output: updated `blast` results with cluster names added; *Figure out naming convention*
	- also picks out representative protein from each cluster `org_prot.bclust.rep.accnums.txt`

### Domain Architectures
#### InterProScan (IPRSCAN)
`iprscan.sb`/`iprscan.sh` <br>
- runs InterProScan analysis on given fasta file(s)
- input: file with fasta sequence(s), takes single fasta sequence or multiple; `org_prot.fa` or `org_prot.all_accnums.fa`
- output: file containing predicted domains w/ standard InterProScan output; `org_prot.iprscan.txt` 

#### ipr2viz
`iprscan_viz.R` <br>
- visualizes iprscan results
- input: iprscan tsv output file
- output: gggenes graph w/ annotated domains

#### ipr2da
`ipr2da.R` <br>
- input: iprscan output, blast results
- output: column with iprscan domains added to blast results
- needs more databases added + commented out

### BLAST+BLASTCLUST+DA

### HEATMAP 
`rep.heatmap.R` <br>
- input: deltablast output
- output: heatmap of homologs to protein of interest in representative organisms (from PATRIC database)
- needs gene name from `xtract_query_info` implemented

### EXTRACTING INFO ABOUT QUERY PROTEIN
`xtract_query_info.R` <br>
- input: accession number of protein of interest


## Output files

### Full analysis

#### Input format
- [x] input.txt with the list of fast file paths

#### Temporary Outputs
- [x] acc2info.tsv
- [x] all_accnums.fa
- [x] all_accnums.txt
- [x] bclust.L60S80.tsv
- [x] blast.cln.tsv
- [x] cln.clust.tsv
- [x] clustIDs
- [x] clust_reps
- [x] full_analysis.tsv
- [x] ipr_domarch.tsv
- [x] iprscan.tsv
- [x] iprscan_cln.tsv
- [ ] rps.out
- [x] wblast.tsv / dblast.tsv

#### Final Outputs (to be read into the app)
- [x] full_analysis.tsv
- [x] iprscan_cln.tsv

### Domain Architecture (DA) analysis

#### Input format
- [x] input.txt with the list of fast file paths

#### Temporary Outputs
- [x] acc2info.tsv
- [x] all_accnums.tsv
- [x] iprscan.tsv

#### Final Outputs (to be read into the app)
- [x] ipr_domarch.tsv
- [x] iprscan_cln.tsv
- output: gene name of protein
- doesn't take any input yet, still working on bugs
