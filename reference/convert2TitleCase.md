# Changing case to 'Title Case'

Translate string to Title Case w/ delimitter.

Translate string to Title Case w/ delimitter. Changing case to 'Title
Case'

## Usage

``` r
convert2TitleCase(text, delimitter)

to_titlecase(text, delimitter)
```

## Arguments

- x:

  Character vector.

- y:

  Delimitter. Default is space (" ").

## Value

Character vector with the input strings converted to title case.

A character vector in title case.

## See also

chartr, toupper, and tolower.

chartr, toupper, and tolower.

## Author

Andrie, Janani Ravi

## Examples

``` r
# Convert a single string to title case
convert2TitleCase("hello world") # Returns "Hello World"
#> [1] "Hello World"

convert2TitleCase("hello world")
#> [1] "Hello World"
convert2TitleCase("this is a test", "_")
#> [1] "This is a test"
```
