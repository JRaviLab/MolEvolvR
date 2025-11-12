# mapAcc2Name

Default renameFA() replacement function. Maps an accession number to its
name

Default rename_fasta() replacement function. Maps an accession number to
its name

## Usage

``` r
mapAcc2Name(line, acc2name, acc_col = "AccNum", name_col = "Name")

mapAcc2Name(line, acc2name, acc_col = "AccNum", name_col = "Name")
```

## Arguments

- line:

  The line of a fasta file starting with '\>'

- acc2name:

  Data Table containing a column of accession numbers and a name column

- acc_col:

  Name of the column containing Accession numbers

- name_col:

  Name of the column containing the names that the accession numbers are
  mapped to

## Value

A character string representing the updated FASTA line, where the
accession number is replaced with its corresponding name.

Character string. The modified line from the Fasta file header with the
name instead of the accession number.

## Examples

``` r
if (FALSE) { # \dontrun{
mapAcc2Name(">P12345 some description", acc2name, "AccNum", "Name")
} # }
if (FALSE) { # \dontrun{
acc2name_table <- data.table(AccNum = c("ACC001", "ACC002"),
Name = c("Species A", "Species B"))
line <- ">ACC001 some additional info"
mapped_line <- mapAcc2Name(line, acc2name_table)
mapped_line  # Expected output: ">Species A"
} # }
```
