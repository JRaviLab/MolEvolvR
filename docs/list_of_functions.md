# List of functions in `the_process` package

## Import and combining input files
- [ ] `scripts/convert_opinscls_tsv.R` using
  - [ ] `clean_clust_file`
  - [ ] add_colnames (currently `colnames.op_ins_cls` and `colnames.op_ins_cls.clus2table`)
  - [ ] remove # rows and convert to columns w/ ID and clust_names
  - [ ] add_uniq_ids (to add GCA_ID, IPG, taxID columns based on AccNums)
  - [ ] add lineage

## Cleanup
- [x] `repeat2s`
  - [x] repeated domains with `(s)`
  - [ ] alternative using `map`?

- [x] `remove_tails`
  - [x] by DomArch
 
- [ ] `remove_empty_rows`
  - [ ] by Species, DomArch, ClustName, GenContext
  
- [ ] `cleanup_clust`
  - [ ] ! `ClustName` acts as `DomArch`?
  - [ ] remove start and end '+'s
  - [ ] `domains_keep` remove rows without query (reads `query_domains`, `domains_keep`)
  - [ ] `domains_rename`   (reads `domains_rename`)
  - [ ] ignored?           (reads `clustnames_ignore`)
  - [ ] `repeat2s`: repeated domains with `(s)`
  - [ ] `remove_tails`
  - [ ] remove empty rows?
  
- [ ] `cleanup_species`
  - [ ] `remove_empty` rows (!change it to an alert about AccNums w/ no lineage/spp)
  - [ ] removes special characters
  - [ ] check if empty rows/taxIDs are because of server retrieval errors!

- [ ] `cleanup_domarch`
  - [ ] ignored           (reads `domains_ignore`)
  - [ ] `domains_keep` remove rows without query (reads `query_domains`, `domains_keep`)
  - [ ] replaced domains  (reads `domains_rename`)
  - [ ] remove start and end '+'s
  - [ ] `repeat2s`: repeated domains with `(s)`
  - [ ] remove empty rows?
  - [ ] `remove_tails`
  
- [ ] `cleanup_gencontext`
  - [x] `reverse_operons`
  - [ ] `repeat2s`: repeated domains with `(s)`
  - [ ] remove empty rows? risky. many eukaryotes don't carry gencontexts!

- [ ] `add_leaves`
  - [ ] `to_titlecase`
  - [ ] `convert_aln2fa` | `convert_aln2tsv`
  - [ ] `add_leaves` | adding leaves based on AccNum, Lineage and Spp.
  - [ ] add DA to `add_leaves` too?
  - [ ] `convert_accnum2fasta`
  - [ ] `filter_for_phylo`

## Summary stats
- [ ] `count_bycol`
- [ ] `generate_wordcount`
  - [ ] `elements2words`
  - [ ] `words2wc`
  - [ ] `filter_freq` | mostly used within other functions
- [ ] `summary_bylin` (for DA, GC)
  - [ ] `summ_DA_byLin`
  - [ ] `summ_GC_byDALin`
  - [ ] `summ_GC_byLin`
- [ ] `summary_stats` (for DA, GC)
  - [ ] `summ_DA`
  - [ ] `summ_GC`
- [x] `total_counts`
- [ ] `find_paralogs`

## Plotting
- [ ] `upset_plot`
  - [ ] should depend on generate_wordcount
- [x] `lineage_DA_plot`
- [x] `lineage_GC_plot`
- [x] `lineage_domain_repeats_plot`??
- [ ] `wordcloud`
- [ ] `msa_tree`
  - [ ] `msa_pdf`
  - [ ] `phylotree`?
- [ ] `prot_network`
  - [ ] by DA/domains
  - [ ] by GC/DA
  - [ ] subsets by prot/domain
