# Contributing to MolEvolvR

This outlines how to propose a change to MolEvolvR.
For a detailed discussion on contributing to this and other tidyverse packages, please see the [development contributing guide](https://rstd.io/tidy-contrib) and our [code review principles](https://code-review.tidyverse.org/).

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file.
This generally means you'll need to edit [roxygen2 comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`, not a `.Rd` file.
You can find the `.R` file that generates the `.Rd` by reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed.
If you’ve found a bug, please file an issue that illustrates the bug with a minimal
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).
See our guide on [how to create a great issue](https://code-review.tidyverse.org/issues/) for more advice.

### Pull request process

- Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis`.

- Install and load the `usethis` package with:

  ```
  install.packages("usethis")

  library("usethis")
  ```

- Clone and fork the MolEvolvR package using:
  ```
  usethis::create_from_github("JRaviLab/MolEvolvR", fork = TRUE)
  ```
- Install BiocManager from Bioconductor:

  ```
  if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(version = "3.19")
  ```

- Install other development dependencies and then ensure that the package passes R CMD check using `devtools`:

  ```
  install.packages("devtools")

  library("devtools")

  devtools::install_dev_deps()

  devtools::check()
  ```

  _If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing._

- Create a Git branch for your pull request (PR). We recommend using:

  ```
  usethis::pr_init("brief-description-of-change")
  ```

- Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
  The title of your PR should briefly describe the change.
  The body of your PR should contain `Fixes #issue-number`.

<!-- *  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>. -->

<!-- No NEWS.md file exists at the moment. Suggest restoring this line when the file exists to avoid a y confusion.-->

### Code style

- New code should follow the tidyverse [style guide](https://style.tidyverse.org).
  You can use the [styler](https://CRAN.R-project.org/package=styler) package to apply these styles, but please don't restyle code that has nothing to do with your PR.
- Lint Your Code: Ensure your code adheres to our style guidelines by using [lintr](https://lintr.r-lib.org/):

  ```
  install.packages("lintr")

  library("lintr")

  lintr::lint("path/to/your/file.R")
  ```

- We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.

- We use [testthat](https://cran.r-project.org/package=testthat) for unit tests.
  Contributions with test cases included are easier to accept.

## Code of Conduct

Please note that the MolEvolvR project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
