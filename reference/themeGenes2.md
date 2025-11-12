# themeGenes2

themeGenes2

## Usage

``` r
themeGenes2()
```

## Value

A ggplot2 theme object.

## Examples

``` r
library(ggplot2)

# Create a sample plot using the custom theme
ggplot(mtcars, aes(x = wt, y = mpg)) +
    geom_point() +
    themeGenes2() +
    labs(title = "Car Weight vs MPG")
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the MolEvolvR package.
#>   Please report the issue at <https://github.com/jravilab/molevolvr/issues>.

```
