
# ConfoundingExplorer

<!-- badges: start -->
<!-- badges: end -->

The `ConfoundingExplorer` package provides a simple shiny app for interactively exploring the effect of confounding between the condition of interest and a batch variable, in terms of the ability to correctly find the variables that are truly differential between the different conditions. The illustrations in the app are performed using simulated data. 

## Installation

You can install `ConfoundingExplorer` from GitHub:

``` r
remotes::install_github("csoneson/ConfoundingExplorer")
```

## Example

To start the app, simply call the `ConfoundingExplorer()` function:

``` r
library(ConfoundingExplorer)
ConfoundingExplorer()
```

