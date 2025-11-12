# plotLineageHeatmap

Generate a lineage plot

## Usage

``` r
plotLineageHeatmap(prot, domains_of_interest, level = 3, label.size = 8)
```

## Arguments

- prot:

  Data frame containing DomArch and Lineage Columns

- domains_of_interest:

  Vector of domains to check for the presence of in all the lineages

- level:

  The max depth of Lineage. ie) i = Kingdom, 2 = Phylum, 3 = class ...

- label.size:

  Size of the text labels

## Value

A ggplot object representing a heatmap (tile plot) of domain repeat
counts across different lineages, with color intensity representing the
occurrence of domains.

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
plotLineageHeatmap(psp_data,
    c(
        "PspA", "Snf7", "Classical-AAA", "PspF", "PspB", "PspC", "ClgR", "PspM",
        "Thioredoxin", "PspN_N", "DUF3046", "LiaI-LiaF-TM", "Toast_rack", "REC",
        "HISKIN", "HAAS", "SHOCT-bihelical", "SHOCT-like", "Tfu_1009", "PspAA",
        "Spermine_synth", "TM-Flotillin", "Band-7", "Betapropeller",
        "MacB_PCD", "FTSW_RODA_SPOVE", "Cest_Tir", "SIGMA-HTH", "GNTR-HTH",
        "DUF2089-HTH", "PadR-HTH", "RHH", "ZnR"
    ),
    level = 2
)
} # }
```
