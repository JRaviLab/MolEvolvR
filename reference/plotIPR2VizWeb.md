# plotIPR2VizWeb

plotIPR2VizWeb

## Usage

``` r
plotIPR2VizWeb(
  infile_ipr,
  accessions,
  analysis = c("Pfam", "Phobius", "TMHMM", "Gene3D"),
  group_by = "Analysis",
  name = "Name",
  text_size = 15,
  legend_name = "ShortName",
  cols = 5,
  rows = 10
)
```

## Arguments

- infile_ipr:

  A path to the input IPR file (TSV format) containing domain
  information.

- accessions:

  A character vector of accession numbers to filter the analysis.

- analysis:

  A character vector specifying the types of analysis to include (e.g.,
  "Pfam", "Phobius", "TMHMM", "Gene3D"). Default is a vector of these
  analyses.

- group_by:

  A string specifying how to group the visualization. Default is
  "Analysis". Options include "Analysis" or "Query".

- name:

  A string representing the name to use for y-axis labels. Default is
  "Name".

- text_size:

  An integer specifying the text size for the plot. Default is 15.

- legend_name:

  A string representing the column to use for legend labels. Default is
  "ShortName".

- cols:

  An integer specifying the number of columns in the facet wrap. Default
  is 5.

- rows:

  An integer specifying the number of rows in the legend. Default is 10.

## Value

A ggplot object representing the domain architecture visualization for
web display.

## Examples

``` r
if (FALSE) { # \dontrun{
plot <- plotIPR2VizWeb(infile_ipr = "path/to/ipr_file.tsv",
                       accessions = c("ACC123", "ACC456"),
                       analysis = c("Pfam", "TMHMM"),
                       group_by = "Analysis",
                       name = "Gene Name",
                       text_size = 15,
                       legend_name = "ShortName",
                       cols = 5,
                       rows = 10)
plot
} # }
```
