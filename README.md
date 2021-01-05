
# ConfoundingExplorer

<!-- badges: start -->
<!-- badges: end -->

The `ConfoundingExplorer` R package provides a simple shiny app for interactively exploring the effect of confounding between a group variable of interest and a batch variable, in terms of the ability to correctly find the variables that are truly differential between the different levels of the group variable. It is mainly intended for teaching purposes and to illustrate important concepts in experimental design, and all analyses are performed using (quite simplistic) simulated data. 

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

The button `Data description` in the lower left corner provides more details about the data generation and the different analysis approaches implemented in the app. 

