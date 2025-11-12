# Lineage Sunburst

Lineage Sunburst

## Usage

``` r
plotLineageSunburst(
  prot,
  lineage_column = "Lineage",
  type = "sunburst",
  levels = 2,
  colors = NULL,
  legendOrder = NULL,
  showLegend = TRUE,
  maxLevels = 5
)
```

## Arguments

- prot:

  Data frame containing a lineage column that the sunburst plot will be
  generated for

- lineage_column:

  String. Name of the lineage column within the data frame. Defaults to
  "Lineage"

- type:

  String, either "sunburst" or "sund2b". If type is "sunburst", a
  sunburst plot of the lineage

- levels:

  Integer. Number of levels the sunburst will have.

- colors:

  A vector of colors for the sunburst plot. If NULL, default colors are
  used.

- legendOrder:

  String vector. The order of the legend. If legendOrder is NULL,

- showLegend:

  Boolean. If TRUE, the legend will be enabled when the component first
  renders.

- maxLevels:

  Integer, the maximum number of levels to display in the sunburst; 5 by
  default, NULL to disable then the legend will be in the descending
  order of the top level hierarchy. will be rendered. If the type is
  sund2b, a sund2b plot will be rendered.

## Value

A sunburst or sund2b plot based on the input lineage data.

## Examples

``` r
if (FALSE) { # \dontrun{
plotLineageSunburst(prot, lineage_column = "Lineage",
type = "sunburst", levels = 3)
} # }
```
