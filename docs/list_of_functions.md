# List of functions in `the_process` package

## Import and combining input files
- [ ] clust2tab
- [ ] filter_clusters
- [ ] clean_clust_file
  - [ ] add_colnames (currently `colnames.op_ins_cls`)
  - [ ] remove # rows and convert to columns w/ ID and clust_names
  - [ ] add_uniq_ids (to add GCA, IPG, taxID columns based on AccNums)

## Cleanup
- [ ] cls_cleanup
  - [ ] remove rows without query (reads `query_names`)
  
- [x] cleanup_species
  - [ ] remove empty rows (change it to an alert about AccNums w/ no lineage/spp)

- [x] replace_doms
  - [ ] ignored           (reads `domains_reomve`)
  - [ ] replaced domains  (reads `domains_rename`)
  - [ ] remove start and end +s
  - [ ] repeated domains with `(s)`

- [x] reverse_operons

- [ ] add_leaves
  - [ ] to_titlecase
  - [ ] adding leaves based on AccNum, Lineage and Spp.
  - [ ] aln2fasta
  - [ ] accnum2fasta
  - [ ] filter_for_phylo

## Summary stats
- [ ] counts_elements
- [ ] elements2words
- [ ] words2wc
- [ ] summary_bylin (for DA, GC)
  - [ ] summ.DA.byLin
  - [ ] summ.GC.byDALin
  - [ ] summ.GC.byLin
- [ ] summary_stats (for DA, GC)
  - [ ] summ.DA
  - [ ] summ.GC
- [x] total_counts
- [ ] find_paralogs

## Plotting
- [x] upset.plot
- [x] lineage.DA.plot
- [x] lineage.neighbors.plot
- [x] lineage.domain_repeats.plot
- [ ] wordcloud
- [ ] msa
- [ ] phylotree
- [ ] msa_tree
- [ ] full network
- [ ] network subsets
