#' Generate data
#'
#' @keywords internal
#'
#' @param nb1c1,nb1c2,nb2c1,nb2c2 Number of samples in the respective groups
#'   (batch 1/2, condition 1/2)
#' @param fraccond,fracbatch Fraction of variables affected by condition and
#'   batch, respectively
#' @param condeffect,batcheffect Condition and batch effect size, respectively
#' @param analysisapproach Type of statistical model to fit
#' @param nvar Number of variables
#' @param seed Random seed
#'
#' @author Charlotte Soneson
#'
#' @return A list with three elements: the simulated data matrix (m), the
#'   results of the statistical test (res) and the sample annotation table
#'   (annot).
#'
#' @importFrom stats rnorm runif lm residuals predict coefficients p.adjust
#'   model.matrix
#' @importFrom limma removeBatchEffect
#'
.generateData <- function(nb1c1, nb1c2, nb2c1, nb2c2, fraccond, fracbatch,
                          condeffect, batcheffect, analysisapproach, nvar,
                          seed) {
    set.seed(seed)
    m <- matrix(stats::rnorm(n = (nb1c1 + nb1c2 + nb2c1 + nb2c2) * nvar,
                             mean = 10, sd = 2),
                nrow = nvar, dimnames = list(
                    paste0("V", seq_len(nvar)),
                    paste0("S", seq_len(nb1c1 + nb1c2 + nb2c1 + nb2c2))))
    cond <- c(rep("C1", nb1c1), rep("C2", nb1c2),
              rep("C1", nb2c1), rep("C2", nb2c2))
    batch <- c(rep("B1", nb1c1), rep("B1", nb1c2),
               rep("B2", nb2c1), rep("B2", nb2c2))
    condvar <- sample(seq_len(nvar), size = fraccond * nvar)
    batchvar <- sample(seq_len(nvar), size = fracbatch * nvar)
    for (i in condvar) {
        m[i, cond == "C2"] <- m[i, cond == "C2"] +
            sign((stats::runif(1) < 0.5) - 0.5) * condeffect
    }
    for (i in batchvar) {
        m[i, batch == "B2"] <- m[i, batch == "B2"] +
            sign((stats::runif(1) < 0.5) - 0.5) * batcheffect
    }
    pvals <- tryCatch({
        if (analysisapproach == "Remove batch effect in advance") {
            m0 <- limma::removeBatchEffect(m, batch = batch)
        } else if (analysisapproach == paste0("Remove batch effect in ",
                                              "advance,\naccounting for ",
                                              "condition")) {
            m0 <- limma::removeBatchEffect(m, batch = batch,
                                           design = stats::model.matrix(~ cond))
        }
        vapply(seq_len(nvar), function(i) {
            if (analysisapproach == "Don't account for batch effect") {
                l <- stats::lm(m[i, ] ~ cond)
            } else if (analysisapproach == "Include batch effect in model") {
                l <- stats::lm(m[i, ] ~ batch + cond)
            } else if (analysisapproach == "Remove batch effect in advance") {
                l <- stats::lm(m0[i, ] ~ cond)
            } else if (analysisapproach == paste0("Remove batch effect in ",
                                                  "advance,\naccounting for ",
                                                  "condition")) {
                l <- stats::lm(m0[i, ] ~ cond)
            } else {
                stop("Not implemented")
            }
            summary(l)$coef["condC2", "Pr(>|t|)"]
        }, NA_real_)},
        error = function(e) rep(NA_real_, nvar))
    res <- data.frame(feature = rownames(m),
                      batchaff = seq_len(nvar) %in% batchvar,
                      condaff = seq_len(nvar) %in% condvar,
                      p.val = pvals,
                      p.adj = stats::p.adjust(pvals, method = "BH"),
                      row.names = rownames(m))
    annot <- data.frame(sample = colnames(m), batch = batch, cond = cond,
                        row.names = colnames(m))
    return(list(m = m, res = res, annot = annot))
}
