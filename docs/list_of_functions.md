# List of functions in `the_process` package

## Import and combining input files
- [ ] clust2tab??
- [ ] filter_clusters??
- [ ] clean_clust_file
  - [ ] add_colnames (currently `colnames.op_ins_cls`)
  - [ ] remove # rows and convert to columns w/ ID and clust_names
  - [ ] add_uniq_ids (to add GCA, IPG, taxID columns based on AccNums)

## Cleanup
- [x] repeat2s
  - [x] repeated domains with `(s)`
  - [ ] alternative using `map`
  
- [ ] cleanup_clust
  - [ ] remove rows without query (reads `query_names`)
  - [ ] repeated domains with `(s)`
  
- [ ] cleanup_species
  - [ ] remove empty rows (change it to an alert about AccNums w/ no lineage/spp)
  - [ ] removes special characters

- [ ] cleanup_domarch
  - [ ] ignored           (reads `domains_reomve`)
  - [ ] replaced domains  (reads `domains_rename`)
  - [ ] remove start and end +s
  - [ ] repeated domains with `(s)`
  - [ ] remove empty rows
  
- [ ] cleanup_gencontext
  - [x] reverse_operons
  - [ ] repeated domains with `(s)`
  - [ ] remove empty rows

- [ ] add_leaves
  - [ ] to_titlecase
  - [ ] adding leaves based on AccNum, Lineage and Spp.
  - [ ] aln2fasta
  - [ ] accnum2fasta
  - [ ] filter_for_phylo

## Summary stats
- [ ] count_bycol
- [ ] generate_wordcount
  - [ ] elements2words
  - [ ] words2wc
- [ ] summary_bylin (for DA, GC)
  - [ ] summ_DA_byLin
  - [ ] summ_GC_byDALin
  - [ ] summ_GC_byLin
- [ ] summary_stats (for DA, GC)
  - [ ] summ_DA
  - [ ] summ_GC
- [x] total_counts
- [ ] find_paralogs

## Plotting
- [ ] upset_plot
  - [ ] should depend on generate_wordcount
- [x] lineage_DA_plot
- [x] lineage_GC_plot
- [x] lineage_domain_repeats_plot??
- [ ] wordcloud
- [ ] msa_tree
  - [ ] msa
  - [ ] phylotree?
- [ ] full network
- [ ] network subsets
