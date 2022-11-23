# Additional functions remaining to be cleaned up
- need to check versions if already in `R`
- if not, need to check usefulness and generalize
- add documentation before moving to `R`

## paralogs
not explored beyond PSP

### compare_paralogs_scripts
- [ ] find_paralogs_new
- [ ] find_paralogs_old
- [ ] original find_paralogs in `R/summarize.R`

## reverter scripts
- not explored beyond PSP

### compare_reverter_scripts
- no functions

### reverse_operons
- `v1` and `v2` -- functions that have not been tested or used beyond PSP

## cleanup
- not explored beyond PSP

### cleanup_domarch
- [ ] no functions but deals with domains_rename, domains_keep, domains_remove

### find_duplicates
- [ ] no functions but removes duplicates to get accurate counts per prot/dom

### identify_da
- [ ] no functions but scripts to identify unique clustnames, domains to retain?

### identify_missing
- [ ] no functions but exploring source of missing data (Species, Lineage, TaxID, GCA_ID, GI, DA, GC) & merging/fixing

### unused_old_cleanup
- [ ] domarch.convert2s.forupset
- [ ] colnames.op_ins_cls
- [ ] repeat2s
- [ ] replace.toastrack
- [ ] remove.empty.rows
- [ ] species.cleanup
- check for newer versions in `R` to figure what to do with these | likely only tested with PSP

## process
- [ ] more recent `Rmd` version in `docs`

## msa_tree

### add_leaves.R
- [ ] to_titlecase | part of `R/pre-msa-tree.R`
- [ ] add_leaves | part of `R/pre-msa-tree.R`

### convert_aln2fa
- [ ] convert_aln2fa | part of `R/pre-msa-tree.R`
- [ ] unused: aln2tsv, convert_to_fasta

### all_aln2fa
- [ ] generate_all_aln2fa | part of `R/pre-msa-tree.R`

### msa-tree
- [ ] no functions but scripts to generate MSA and phylogenetic tree
- [ ] check to see if it overlaps with `R/pre-msa-tree`

### tree
- [ ] seq_tree | file --> tree
- [ ] few other functions that don't work!
- newer/modified version in `R/tree.R`? check!

## summarize/dataviz
### figx-xxx-analysis-plots
- [ ] no functions but ... similar to `the_process`
  - [ ] data cleanup | Species, Lineage, DomArch, GenContext
  - [ ] counts | DA, GC
  - [ ] word counts, wordcloud
  - [ ] upSet | da-doms
  - [ ] stacked lineage plots

### upsetr
- [ ] with directionality | radlinsky | unsused? check!

### network
- `20190408`, `av`, `ngrams` -- no functions, and maybe old and unused

### total_counts
- quick script to count uniq GC (only tested on PSP)

## Missing Functions
Including the specific use case ones from e.g., PSP <br>
- `find_duplicates`
- `fill_missing`
