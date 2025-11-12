# InterProScan Column Names

A character vector containing the expected column names from an
InterProScan output table. This dataset is useful for validating,
parsing, or reconstructing data frames produced by InterProScan.

## Usage

``` r
ipr_colnames
```

## Format

A character vector with 13 elements:

- AccNum:

  Accession number of the sequence.

- SeqMD5Digest:

  MD5 digest of the sequence.

- SLength:

  Length of the sequence.

- Analysis:

  Type of analysis or database used (e.g., Pfam, SMART).

- DB.ID:

  Database-specific identifier.

- SignDesc:

  Description of the signature or domain.

- StartLoc:

  Start position of the match on the sequence.

- StopLoc:

  Stop position of the match on the sequence.

- Score:

  Score assigned to the match (if applicable).

- Status:

  Status of the analysis (e.g., OK, WARNING).

- RunDate:

  Date the InterProScan analysis was run.

- IPRAcc:

  InterPro accession number.

- IPRDesc:

  InterPro entry description.

## Source

Generated internally to represent standard InterProScan output fields.

## Examples

``` r
data(ipr_colnames)
ipr_colnames
#>  [1] "AccNum"       "SeqMD5Digest" "SLength"      "Analysis"     "DB.ID"       
#>  [6] "SignDesc"     "StartLoc"     "StopLoc"      "Score"        "Status"      
#> [11] "RunDate"      "IPRAcc"       "IPRDesc"     
```
