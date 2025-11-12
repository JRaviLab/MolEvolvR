# addLineage

addLineage

addLineage

## Usage

``` r
addLineage(
  df,
  acc_col = "AccNum",
  assembly_path,
  lineagelookup_path,
  ipgout_path = NULL,
  plan = "multicore"
)

addLineage(
  df,
  acc_col = "AccNum",
  assembly_path,
  lineagelookup_path,
  ipgout_path = NULL,
  plan = "multicore"
)
```

## Arguments

- df:

  Dataframe containing accession numbers. The dataframe should have a
  column specified by `acc_col` that contains these accession numbers.

- acc_col:

  Character. The name of the column in `df` containing accession
  numbers. Default is "AccNum".

- assembly_path:

  String. The path to the assembly summary file generated using the
  [`downloadAssemblySummary()`](https://jravilab.github.io/MolEvolvR/reference/downloadAssemblySummary.md)
  function.

- lineagelookup_path:

  String. The path to the lineage lookup file (taxid to lineage mapping)
  generated using the `create_lineage_lookup()` function.

- ipgout_path:

  String. Optional path to save intermediate output files. Default is
  NULL.

- plan:

  Character. Specifies the execution plan for parallel processing.
  Default is "multicore".

## Value

A `data.frame` that combines the original `df` with the lineage
information.

A dataframe that combines the original dataframe `df` with lineage
information retrieved based on the provided accession numbers.

## Examples

``` r
if (FALSE) { # \dontrun{
addLineage()
} # }
if (FALSE) { # \dontrun{
enriched_df <- addLineage(df = my_data,
                           acc_col = "AccNum",
                           assembly_path = "path/to/assembly_summary.txt",
                           lineagelookup_path = "path/to/lineage_lookup.tsv")
} # }
```
