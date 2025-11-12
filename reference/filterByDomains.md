# filterByDomains

filterByDomains filters a data frame by identifying exact domain matches
and either keeping or removing rows with the identified domain

## Usage

``` r
filterByDomains(
  prot,
  column = "DomArch",
  doms_keep = c(),
  doms_remove = c(),
  ignore.case = FALSE
)
```

## Arguments

- prot:

  Dataframe to filter

- column:

  Column to search for domains in (DomArch column)

- doms_keep:

  Vector of domains that must be identified within column in order for
  observation to be kept

- doms_remove:

  Vector of domains that, if found within an observation, will be
  removed

- ignore.case:

  Should the matching be non case sensitive

## Value

Filtered data frame

## Note

There is no need to make the domains 'regex safe', that will be handled
by this function

## Author

Samuel Chen, Janani Ravi

## Examples

``` r
if (FALSE) { # \dontrun{
filterByDomains()
} # }
```
