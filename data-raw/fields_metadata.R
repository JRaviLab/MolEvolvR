## code to prepare `fields_metadata` dataset goes here

fields_metadata <- list(
  submission_type=list(
    name="Submission Type",
    values=list(
      'full'="FASTA",
      'phylo'="Phylogenetic Analysis",
      'da'="Domain Architecture",
      'blast'="BLAST",
      'dblast'="BLAST"
    )
  ),
  homology_search=list(
    name="Homology Search?"
  ),
  submitter_email=list(
    name="Submitter Email"
  ),
  database=list(
    name="Database for homology search"
  ),
  nhits=list(
    name="Maximum hits"
  ),
  evalue=list(
    name="E-value cutoff"
  ),
  includes_ncbi_acc=list(
    name="Contains NCBI accessions?",
    values=list(
      "TRUE"="Yes",
      "FALSE"="No"
    )
  ),
  advanced_options=list(
    name="Advanced Options",
    values=list(
      "domain_architecture"="Domain Architecture",
      "homology_search"="Homology Search",
      "phylogenetic_analysis"="Phylogenetic Analysis"
    )
  ),
  job_code=list(
    name="Job Code"
  )
)

usethis::use_data(fields_metadata, overwrite = TRUE)
