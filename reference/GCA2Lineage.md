# Function to map GCA_ID to TaxID, and TaxID to Lineage

Function to map GCA_ID to TaxID, and TaxID to Lineage

## Usage

``` r
GCA2Lineage(
  prot_data,
  assembly_path = "/data/research/jravilab/common_data/assembly_summary_genbank.txt",
  lineagelookup_path = "/data/research/jravilab/common_data/lineage_lookup.tsv",
  acc_col = "AccNum"
)
```

## Arguments

- prot_data:

  Dataframe containing a column `GCA_ID`

- assembly_path:

  String of the path to the assembly_summary path This file can be
  generated using the "downloadAssemblySummary()" function

- lineagelookup_path:

  String of the path to the lineage lookup file (taxid to lineage
  mapping). This file can be generated using the "createLineageLookup()"
  function

- acc_col:

  Character. The name of the column in `prot_data` containing accession
  numbers. Default is "AccNum".

## Value

A dataframe containing the merged information of GCA_IDs, TaxIDs, and
their corresponding lineage up to the phylum level. The dataframe will
include information from the input `prot_data` and lineage data.

## Note

Currently configured to have at most kingdom and phylum

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
result <- GCA2Lineage(prot_data = my_prot_data,
                       assembly_path = "path/to/assembly_summary.txt",
                       lineagelookup_path = "path/to/lineage_lookup.tsv")
} # }
```
