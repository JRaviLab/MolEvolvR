## Test script to look up assembly accession based on protein accession numbers
## using IPG database and edirect, efetch

## Loading modules
module purge
module use /mnt/home/johnj/software/modulefiles
module load R/4.0.2
module load edirect

# efetch + prot accnum --> ipg
printf "\ntesting single accnum\n"
efetch -db protein -id "CBE06962.1" -format ipg
## with a list
printf "\n\ntesting list of accnum\n"
## print statement for current step
for x in acclist.txt                                       
do
   accnum=$(cat $x | tr '\n' ',') # since efetch likes commas instead of \n char                           
   efetch -db protein -format ipg -id "$accnum" >> test_prot2ipg.out
   # this output file can then be stitched back to the original BLAST one
done

