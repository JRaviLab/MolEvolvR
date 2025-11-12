# plotWordCloud3

plotWordCloud3

## Usage

``` r
wordcloud3(
  data,
  size = 1,
  minSize = 0,
  gridSize = 0,
  fontFamily = "Segoe UI",
  fontWeight = "bold",
  color = "random-dark",
  backgroundColor = "white",
  minRotation = -pi/4,
  maxRotation = pi/4,
  shuffle = TRUE,
  rotateRatio = 0.4,
  shape = "circle",
  ellipticity = 0.65,
  widgetsize = NULL,
  figPath = NULL,
  hoverFunction = NULL
)
```

## Arguments

- data:

  Data frame or table containing words and their frequencies for the
  word cloud.

- size:

  Numeric. Scaling factor for word sizes (default is 1).

- minSize:

  Numeric. Minimum font size for the smallest word (default is 0).

- gridSize:

  Numeric. Size of the grid for placing words (default is 0).

- fontFamily:

  Character. Font family to use for the words (default is "Segoe UI").

- fontWeight:

  Character. Font weight for the words (default is "bold").

- color:

  Character or vector. Color of the words. Use "random-dark" for random
  dark colors (default) or specify a color.

- backgroundColor:

  Character. Background color of the word cloud (default is "white").

- minRotation:

  Numeric. Minimum rotation angle of words in radians (default is -π/4).

- maxRotation:

  Numeric. Maximum rotation angle of words in radians (default is π/4).

- shuffle:

  Logical. Whether to shuffle the words (default is TRUE).

- rotateRatio:

  Numeric. Proportion of words that are rotated (default is 0.4).

- shape:

  Character. Shape of the word cloud ("circle" is default, but you can
  use "cardioid", "star", "triangle", etc.).

- ellipticity:

  Numeric. Degree of ellipticity (default is 0.65).

- widgetsize:

  Numeric vector. Width and height of the widget (default is NULL, which
  uses default size).

- figPath:

  Character. Path to an image file to use as a mask for the word cloud
  (optional).

- hoverFunction:

  JS function. JavaScript function to run when hovering over words
  (optional).

## Value

An HTML widget object displaying a word cloud.

## Examples

``` r
if (FALSE) { # \dontrun{
wordcloud3(data = your_data, size = 1.5, color = "random-light")
} # }
```
