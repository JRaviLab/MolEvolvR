# Download the combined assembly summaries of genbank and refseq

Download the combined assembly summaries of genbank and refseq

## Usage

``` r
downloadAssemblySummary(
  outpath,
  keep = c("assembly_accession", "taxid", "species_taxid", "organism_name")
)
```

## Arguments

- outpath:

  String of path where the assembly summary file should be written

- keep:

  Character vector containing which columns should be retained and
  downloaded

## Value

A tab-separated file containing the assembly summary. The function does
notreturn any value but writes the output directly to the specified
file.

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
downloadAssemblySummary(outpath = "assembly_summary.tsv",
     keep = c("assembly_accession", "taxid", "organism_name"))
} # }
```
