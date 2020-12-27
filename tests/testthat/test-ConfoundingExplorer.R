test_that("an app is generated", {
  res <- ConfoundingExplorer()
  expect_is(res, "shiny.appobj")
})
