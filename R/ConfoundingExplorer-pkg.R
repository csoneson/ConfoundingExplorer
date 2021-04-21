#' ConfoundingExplorer
#'
#' ConfoundingExplorer is an R package for exploring the impact of confounding
#' between a condition of interest and a batch variable.
#'
#' @name ConfoundingExplorer-pkg
#' @docType package
NULL

globalVariables(c("X", "Y", "jitteredX", "jitteredY",
                  "p.val", "p.adj", "XWidth", "YWidth",
                  "batch", "batchaff", "cond", "condaff", "feature", "value"))
