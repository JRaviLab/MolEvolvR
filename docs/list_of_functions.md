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
 
- [x] `remove_empty_rows`
  - [x] by Species, DomArch, ClustName, GenContext
  
- [x] `cleanup_clust`
  - [ ] ! `ClustName` acts as `DomArch`?
  - [x] remove start and end '+'s
  - [x] `domains_keep` remove rows without query (reads `query_domains`, `domains_keep`)
  - [x] `domains_rename`   (reads `domains_rename`)
  - [ ] ignored?           (reads `clustnames_ignore`)
  - [x] `repeat2s`: repeated domains with `(s)`
  - [x] `remove_tails`
  - [x] remove empty rows?
  
- [x] `cleanup_species`
  - [ ] `remove_empty` rows (!change it to an alert about AccNums w/ no lineage/spp)
  - [x] removes special characters
  - [x] check if empty rows/taxIDs are because of server retrieval errors!

- [x] `cleanup_domarch`
  - [ ] ignored           (reads `domains_ignore`)
  - [x] `domains_keep` remove rows without query (reads `query_domains`, `domains_keep`)
  - [x] replaced domains  (reads `domains_rename`)
  - [x] remove start and end '+'s
  - [x] `repeat2s`: repeated domains with `(s)`
  - [x] remove empty rows?
  - [x] `remove_tails`
  
- [x] `cleanup_gencontext`
  - [x] `reverse_operons`
  - [x] `repeat2s`: repeated domains with `(s)`
  - [ ] remove empty rows? risky. many eukaryotes don't carry gencontexts!

- [x] `add_leaves`
  - [x] `to_titlecase`
  - [x] `convert_aln2fa` | `convert_aln2tsv`
  - [x] `add_leaves` | adding leaves based on AccNum, Lineage and Spp.
  - [ ] add DA to `add_leaves` too?
  - [ ] `convert_accnum2fasta`
  - [ ] `filter_for_phylo`

## Summary stats
- [x] `count_bycol`
- [ ] `generate_wordcount`
  - [x] `elements2words`
  - [x] `words2wc`
  - [ ] `filter_freq` | mostly used within other functions
- [x] `summary_bylin` (for DA, GC)
  - [x] `summ_DA_byLin`
  - [x] `summ_GC_byDALin`
  - [x] `summ_GC_byLin`
- [x] `summary_stats` (for DA, GC)
  - [x] `summ_DA`
  - [x] `summ_GC`
- [x] `total_counts`
- [x] `find_paralogs`

## Plotting
- [x] `upset_plot`
  - [ ] should depend on generate_wordcount
- [x] `lineage_DA_plot`
- [x] `lineage_GC_plot`
- [x] `lineage_domain_repeats_plot`??
- [x] `wordcloud`
- [ ] `msa_tree`
  - [ ] `msa_pdf`
  - [ ] `phylotree`?
- [x] `prot_network`
  - [x] by DA/domains
  - [x] by GC/DA
  - [x] subsets by prot/domain
