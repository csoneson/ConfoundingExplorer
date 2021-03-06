% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConfoundingExplorer.R
\name{ConfoundingExplorer}
\alias{ConfoundingExplorer}
\title{Confounding explorer}
\usage{
ConfoundingExplorer(
  sampleSizes = matrix(rep(5, 4), nrow = 2, dimnames = list(c("group1", "group2"),
    c("batch1", "batch2"))),
  fracVarCond = 0.25,
  fracVarBatch = 0.5,
  fracVarUnknown = 0,
  condEffectSize = 3,
  batchEffectSize = 3,
  unknownEffectSize = 0,
  unknownEffectType = "categorical",
  analysisApproach = "dontAdjust",
  seed = 123
)
}
\arguments{
\item{sampleSizes}{2x2 numeric matrix giving the number of samples in each
group. Row names must be c('group1', 'group2') and column names must be
c('batch1', 'batch2').}

\item{fracVarCond, fracVarBatch, fracVarUnknown}{Numeric scalars between 0 and
1. The fraction of variables affected by the condition effect, batch
effect, and 'unknown' effect, respectively.}

\item{condEffectSize, batchEffectSize, unknownEffectSize}{Numeric scalars.
The condition, batch and 'unknown' effect size, respectively.}

\item{unknownEffectType}{Character scalar, either 'categorical' or
'continuous', representing the type of 'unknown' effect to add.}

\item{analysisApproach}{Character scalar. One of 'dontAdjust', 'inclBatch',
'removeBatch', 'removeBatchAccCond'. Determines what model is fit to the
data.}

\item{seed}{Numeric scalar, the random seed to use when simulating data.}
}
\value{
A shinyApp object
}
\description{
Confounding explorer
}
\examples{
if (interactive()) {
  ConfoundingExplorer()
}

}
\author{
Charlotte Soneson
}
