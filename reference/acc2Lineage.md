# acc2Lineage

This function combines 'efetchIPG()' and 'IPG2Lineage()' to map a set of
protein accessions to their assembly (GCA_ID), tax ID, and lineage.

Function to map protein accession numbers to lineage

This function combines 'efetchIPG()' and 'IPG2Lineage()' to map a set of
protein accessions to their assembly (GCA_ID), tax ID, and lineage.

## Usage

``` r
acc2Lineage(
  accessions,
  assembly_path,
  lineagelookup_path,
  ipgout_path = NULL,
  plan = "multicore"
)

acc2Lineage(
  accessions,
  assembly_path,
  lineagelookup_path,
  ipgout_path = NULL,
  plan = "multicore"
)
```

## Arguments

- accessions:

  Character vector of protein accessions

- assembly_path:

  String of the path to the assembly_summary path This file can be
  generated using the "downloadAssemblySummary()" function

- lineagelookup_path:

  String of the path to the lineage lookup file (taxid to lineage
  mapping). This file can be generated using the

- ipgout_path:

  Path to write the results of the efetch run of the accessions on the
  ipg database. If NULL, the file will not be written. Defaults to NULL

- plan:

  Character. Specifies the execution plan for parallel processing.
  Default is "multicore".

## Value

A `data.table` that contains the lineage information, mapping protein
accessions to their tax IDs and lineages.

A dataframe containing lineage information mapped to the given protein
accessions. The dataframe includes relevant columns such as TaxID,
GCA_ID, Protein, Protein Name, Species, and Lineage.

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
acc2Lineage()
} # }
if (FALSE) { # \dontrun{
lineage_data <- acc2Lineage(
  accessions = c("P12345", "Q67890"),
  assembly_path = "path/to/assembly_summary.txt",
  lineagelookup_path = "path/to/lineage_lookup.tsv",
  ipgout_path = "path/to/output.txt"
)
} # }
```
