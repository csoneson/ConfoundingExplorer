test_that("an app is generated", {
  res <- ConfoundingExplorer()
  expect_s3_class(res, "shiny.appobj")
})

test_that("data generation works", {
    analysis_approaches <- c(
        "Don't account for batch effect",
        "Include batch effect in model",
        "Remove batch effect in advance",
        "Remove batch effect in advance,\naccounting for condition"
    )

    for (appr in analysis_approaches) {
        dat <- .generateData(nb1c1 = 5, nb1c2 = 5, nb2c1 = 5, nb2c2 = 5,
                             fraccond = 0.5, fracbatch = 0.5,
                             condeffect = 2, batcheffect = 2,
                             analysisapproach = appr,
                             nvar = 100, seed = 123)
        expect_type(dat, "list")
        expect_named(dat, c("m", "res", "annot"))
        expect_equal(sum(dat$res$batchaff), 0.5 * 100)
        expect_equal(sum(dat$res$condaff), 0.5 * 100)
    }
})

test_that("no p-values are calculated if adjusting with full confounding", {
    dat <- .generateData(nb1c1 = 5, nb1c2 = 0, nb2c1 = 0, nb2c2 = 5,
                         fraccond = 0.5, fracbatch = 0.5,
                         condeffect = 2, batcheffect = 2,
                         analysisapproach = "Include batch effect in model",
                         nvar = 100, seed = 123)
    expect_type(dat, "list")
    expect_named(dat, c("m", "res", "annot"))
    expect_equal(dat$res$p.val, rep(NA_real_, length(dat$res$p.val)))
    expect_equal(dat$res$p.adj, rep(NA_real_, length(dat$res$p.adj)))
})
