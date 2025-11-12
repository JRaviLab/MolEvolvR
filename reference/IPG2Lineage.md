# IPG2Lineage

Takes the resulting file of an efetch run on the ipg database and

Takes the resulting file of an efetch run on the ipg database and append
lineage, and taxid columns

## Usage

``` r
IPG2Lineage(
  accessions,
  ipg_file,
  refseq_assembly_path,
  genbank_assembly_path,
  lineagelookup_path
)

IPG2Lineage(
  accessions,
  ipg_file,
  refseq_assembly_path,
  genbank_assembly_path,
  lineagelookup_path
)
```

## Arguments

- accessions:

  Character vector of protein accessions

- ipg_file:

  Path to the file containing results of an efetch run on the ipg
  database. The protein accession in 'accessions' should be contained in
  this file

- refseq_assembly_path:

  String. Path to the RefSeq assembly summary file.

- genbank_assembly_path:

  String. Path to the GenBank assembly summary file.

- lineagelookup_path:

  String of the path to the lineage lookup file (taxid to lineage
  mapping). This file can be generated using the "createLineageLookup()"
  function

- assembly_path:

  String of the path to the assembly_summary path This file can be
  generated using the
  [downloadAssemblySummary](https://jravilab.github.io/MolEvolvR/reference/downloadAssemblySummary.md)
  function

## Value

A `data.table` with the lineage information for the provided protein
accessions.

A data table containing protein accessions along with their
corresponding TaxIDs and lineage information.

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
IPG2Lineage()
} # }

if (FALSE) { # \dontrun{
lins <- IPG2Lineage(
  accessions = c("P12345", "Q67890"),
  ipg_file = "path/to/ipg_results.txt",
  refseq_assembly_path = "path/to/refseq_assembly_summary.txt",
  genbank_assembly_path = "path/to/genbank_assembly_summary.txt",
  lineagelookup_path = "path/to/lineage_lookup.tsv"
)
} # }
```
