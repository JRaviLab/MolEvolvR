source("./R/pre-msa-tree.R")
source("./molevol_scripts/filter_tree.R")
source("scripts/tree.R")

acc = filter_tree("../molevol_data/project_data/phage_defense/full_analysis_20210108/cln_combined.tsv",
                  domains_of_interest =  c("P-loop_containing_nucleotide_triphosphate_hydrolases",
                                           "Cytidine_Deaminase_domain_2"),
                  interest_col = "DomArch.Gene3D", tail_cutoff = 5)

cln_combined = fread("../molevol_data/project_data/phage_defense/full_analysis_20210108/cln_combined.tsv",
                     sep = "\t", fill = T)

tmp_fa = tempfile()
acc2fa(acc, tmp_fa)

### Rename
renamed = rename_fasta(fa_path = tmp_fa,
                       outpath = tmp_fa, replacement_function = map_acc2name,
                       acc2name = cln_combined)

view(renamed)

tmp_msa = tempfile()

alignFasta(tmp_fa, tool = "ClustalOmega" ,outpath = tmp_msa )

seq_tree(tmp_msa)
