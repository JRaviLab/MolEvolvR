# MolEvolvR Summary

A collection of summary functions for the MolEvolvR package.

Function to summarize and retrieve counts by Domains & Domains+Lineage

Function to retrieve counts of how many lineages a DomArch appears in

Creates a data frame with a totalcount column

This function is designed to sum the counts column by either Genomic
Context or Domain Architecture and creates a totalcount column from
those sums.

## Usage

``` r
summarizeByLineage(prot = "prot", column = "DomArch", by = "Lineage", query)

summarizeDomArch_ByLineage(x)

summarizeDomArch(x)

summarizeGenContext_ByDomArchLineage(x)

summarizeGenContext_ByLineage(x)

summarizeGenContext(x)

totalGenContextOrDomArchCounts(
  prot,
  column = "DomArch",
  lineage_col = "Lineage",
  cutoff = 90,
  RowsCutoff = FALSE,
  digits = 2
)
```

## Arguments

- prot:

  A data frame that must contain columns:

  - count

- column:

  Character. The column to summarize, default is "DomArch".

- by:

  A string representing the grouping column (e.g., `Lineage`). Default
  is "Lineage".

- query:

  A string specifying the query pattern for filtering the target column.
  Use "all" to skip filtering and include all rows.

- x:

  A dataframe or tibble containing the data. It must have columns named
  `GenContext`, `DomArch`, `Lineage`, and `count`.

- lineage_col:

  Character. The name of the lineage column, default is "Lineage".

- cutoff:

  Numeric. Cutoff for total count. Counts below this cutoff value will
  not be shown. Default is 0.

- RowsCutoff:

  Logical. If TRUE, filters based on cumulative percentage cutoff.
  Default is FALSE.

- digits:

  Numeric. Number of decimal places for percentage columns. Default is
  2.

## Value

A tibble summarizing the counts of occurrences of elements in the
`column`, grouped by the `by` column. The result includes the number of
occurrences (`count`) and is arranged in descending order of count.

A tibble summarizing the counts of unique domain architectures
(`DomArch`) per lineage (`Lineage`). The resulting table contains three
columns: `DomArch`, `Lineage`, and `count`, which indicates the
frequency of each domain architecture for each lineage. The results are
arranged in descending order of `count`.

A tibble summarizing each unique `DomArch`, along with the following
columns:

- `totalcount`: The total occurrences of each `DomArch` across all
  lineages.

- `totallin`: The total number of unique lineages in which each
  `DomArch` appears. The results are arranged in descending order of
  `totallin` and `totalcount`.

A tibble summarizing each unique combination of `GenContext`, `DomArch`,
and `Lineage`, along with the following columns:

- `GenContext`: The genomic context for each entry.

- `DomArch`: The domain architecture for each entry.

- `Lineage`: The lineage associated with each entry.

- `count`: The total number of occurrences for each combination of
  `GenContext`, `DomArch`, and `Lineage`.

The results are arranged in descending order of `count`.

A tibble summarizing each unique combination of `GenContext` and
`Lineage`, along with the count of occurrences. The results are arranged
in descending order of count.

A tibble summarizing each unique `GenContext`, along with the following
columns:

- `totalcount`: The total count for each `GenContext`.

- `totalDA`: The number of distinct `DomArch` for each `GenContext`.

- `totallin`: The number of distinct `Lineage` for each `GenContext`.

The results are arranged in descending order of `totalcount`, `totalDA`,
and `totallin`.

A data frame with the following columns:

- `{{ column }}`: Unique values from the specified column.

- `totalcount`: The total count of occurrences for each unique value in
  the specified column.

- `IndividualCountPercent`: The percentage of each `totalcount` relative
  to the overall count.

- `CumulativePercent`: The cumulative percentage of total counts.

## Note

Please refer to the source code if you have alternate file formats
and/or column names.

## Examples

``` r
if (FALSE) { # \dontrun{
library(tidyverse)
tibble(DomArch = c("a+b", "a+b", "b+c", "a+b"), Lineage = c("l1", "l1", "l1", "l2")) |>
    summarizeByLineage(query = "all")
} # }

if (FALSE) { # \dontrun{
summarizeDomArch_ByLineage(data1)
} # }
if (FALSE) { # \dontrun{
summarizeDomArch(data1)
} # }
if (FALSE) { # \dontrun{
summarizeGenContext_ByDomArchLineage(your_data)
} # }
if (FALSE) { # \dontrun{
summarizeGenContext_ByLineage(your_data)
} # }
if (FALSE) { # \dontrun{
summarizeGenContext(data1)
} # }
if (FALSE) { # \dontrun{
totalGenContextOrDomArchCounts(pspa - gc_lin_counts, 0, "GC")
} # }
```
