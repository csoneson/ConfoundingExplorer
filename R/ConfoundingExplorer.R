#' Confounding explorer
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @importFrom shiny sliderInput radioButtons numericInput fluidRow column
#'   reactiveValues observe validate need renderPlot shinyApp
#' @importFrom ggplot2 ggplot geom_tile aes geom_point labs scale_x_discrete
#'   scale_y_discrete scale_color_manual theme_bw theme element_text
#'   geom_histogram
#' @importFrom iSEE jitterSquarePoints
#' @importFrom iCOBRA COBRAData calculate_adjp calculate_performance
#'   prepare_data_for_plot plot_fpr plot_tpr
#' @importFrom cowplot plot_grid
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#'   menuItem dashboardBody box
#' @importFrom ComplexHeatmap columnAnnotation rowAnnotation Heatmap
#' @importFrom shinyMatrix matrixInput
#' @importFrom grDevices rgb
#' @importFrom stats coefficients lm p.adjust predict residuals rnorm runif
#'
#' @return A shinyApp object
#'
#' @examples
#' if (interactive) {
#'   ConfoundingExplorer()
#' }
ConfoundingExplorer <- function() {
    ui <- shinydashboard::dashboardPage(
        skin = "red",

        # Application title
        header = shinydashboard::dashboardHeader(
            title = "Experimental design and confounding exploration",
            titleWidth = 600
        ),

        # Sidebar
        sidebar = shinydashboard::dashboardSidebar(
            width = 300,

            shinydashboard::menuItem(
                "Set group sizes",
                startExpanded = TRUE,

                shinyMatrix::matrixInput(
                    "sampleSizeTable",
                    class = "numeric",
                    label = "Group sizes",
                    rows = list(
                        names = TRUE,
                        editableNames = FALSE),
                    cols = list(
                        names = TRUE,
                        editableNames = FALSE),
                    value = matrix(rep(5, 4), nrow = 2,
                                   dimnames = list(c("group1", "group2"),
                                                   c("batch1", "batch2")))
                )

            ),
            shinydashboard::menuItem(
                "Set fraction of affected variables",
                startExpanded = TRUE,

                shiny::sliderInput(
                    "fraccond",
                    "Fraction of variables affected by condition:",
                    min = 0, max = 1, value = 0.25
                ),
                shiny::sliderInput(
                    "fracbatch",
                    "Fraction of variables affected by batch:",
                    min = 0, max = 1, value = 0.5
                )
            ),
            shinydashboard::menuItem(
                "Set effect sizes",
                startExpanded = TRUE,

                shiny::sliderInput(
                    "condeffect", "Condition effect size:",
                    min = 0, max = 10, step = 0.1, value = 3
                ),
                shiny::sliderInput(
                    "batcheffect", "Batch effect size:",
                    min = 0, max = 10, step = 0.1, value = 3
                )
            ),
            shiny::radioButtons(
                "analysisapproach",
                "Analysis approach",
                choices = c("Don't account for batch effect",
                            "Include batch effect in model",
                            "Remove batch effect in advance",
                            paste0("Remove batch effect in advance,",
                                   "\n, accounting for condition")),
                selected = "Include batch effect in model"
            ),

            shiny::numericInput(
                "seed", "Random seed", 123, min = 1, max = 1e8
            ),

            shiny::actionButton("datadesc", "Data description")
        ),

        body = shinydashboard::dashboardBody(
            shiny::fluidRow(
                shiny::column(
                    width = 6,
                    shinydashboard::box(
                        title = "Statistical test results (Condition effect)",
                        width = NULL,
                        shiny::plotOutput("hintonPlot")
                    )
                ),
                shiny::column(
                    width = 6,
                    shinydashboard::box(
                        title = "P-value histogram (Condition effect)",
                        width = NULL,
                        shiny::plotOutput("pvalHist")
                    )
                )
            ),
            shiny::fluidRow(
                shiny::column(
                    width = 12,
                    shinydashboard::box(
                        title = "Data heatmap",
                        width = NULL,
                        shiny::plotOutput("heatmapPlot")
                    )
                )
            )
        )
    )

    generateData <- function(nb1c1, nb1c2, nb2c1, nb2c2, fraccond, fracbatch,
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
            vapply(seq_len(nvar), function(i) {
                if (analysisapproach == "Don't account for batch effect") {
                    l <- stats::lm(m[i, ] ~ cond)
                } else if (analysisapproach == "Include batch effect in model") {
                    l <- stats::lm(m[i, ] ~ batch + cond)
                } else if (analysisapproach == "Remove batch effect in advance") {
                    l0 <- stats::lm(m[i, ] ~ batch)
                    m0 <- stats::residuals(l0) + stats::coefficients(l0)["(Intercept)"]
                    l <- stats::lm(m0 ~ cond)
                } else if (analysisapproach == "Remove batch effect in advance,\n, accounting for condition") {
                    l0 <- stats::lm(m[i, ] ~ cond)
                    l1 <- stats::lm(stats::residuals(l0) ~ batch)
                    m0 <- stats::predict(l0) + stats::residuals(l1) +
                        stats::coefficients(l0)["(Intercept)"]
                    l <- stats::lm(m0 ~ cond)
                } else {
                    stop("Not implemented")
                }
                summary(l)$coef["condC2", "Pr(>|t|)"]
            }, NA_real_)},
            error = function(e) rep(NA, nvar))
        res <- data.frame(feature = rownames(m),
                          batchaff = seq_len(nvar) %in% batchvar,
                          condaff = seq_len(nvar) %in% condvar,
                          pval = pvals,
                          padj = stats::p.adjust(pvals, method = "BH"),
                          row.names = rownames(m))
        annot <- data.frame(sample = colnames(m), batch = batch, cond = cond,
                            row.names = colnames(m))
        return(list(m = m, res = res, annot = annot))
    }

    server <- function(input, output, session) {

        ## Display data generation description
        shiny::observeEvent(input$datadesc, {
            shiny::showModal(shiny::modalDialog(
                title = "Data generation",
                shiny::renderUI(shiny::HTML(readLines(
                    system.file("extdata/dataDescription.html",
                                package = "ConfoundingExplorer")))),
                # shiny::includeMarkdown(),
                easyClose = TRUE
            ))
        })

        ## Generate data
        datres <- shiny::reactiveValues()
        shiny::observe({
            shiny::validate(
                shiny::need(is.numeric(input$seed),
                            "Seed must be a number")
            )
            tmp <- generateData(nb1c1 = input$sampleSizeTable[1, 1],
                                nb1c2 = input$sampleSizeTable[2, 1],
                                nb2c1 = input$sampleSizeTable[1, 2],
                                nb2c2 = input$sampleSizeTable[2, 2],
                                fraccond = input$fraccond,
                                fracbatch = input$fracbatch,
                                condeffect = input$condeffect,
                                batcheffect = input$batcheffect,
                                analysisapproach = input$analysisapproach,
                                nvar = 1000, seed = input$seed)
            datres$res <- tmp$res
            datres$m <- tmp$m
            datres$annot <- tmp$annot
        })

        output$heatmapPlot <- shiny::renderPlot({
            shiny::validate(
                shiny::need(!is.null(datres$m) & !is.null(datres$res),
                            "Results could not be generated")
            )
            colAnnot <- ComplexHeatmap::columnAnnotation(
                batch = datres$res$batchaff,
                cond = datres$res$condaff,
                col = list(batch = c(`TRUE` = "forestgreen", `FALSE` = "grey85"),
                           cond = c(`TRUE` = "purple", `FALSE` = "grey85"))
            )
            rowAnnot <- ComplexHeatmap::rowAnnotation(
                batch = datres$annot$batch,
                cond = datres$annot$cond,
                col = list(
                    batch = c(
                        B1 = grDevices::rgb(181, 24, 25, maxColorValue = 256),
                        B2 = grDevices::rgb(37, 202, 67, maxColorValue = 256)
                    ),
                    cond = c(
                        C1 = grDevices::rgb(216, 44, 168, maxColorValue = 256),
                        C2 = grDevices::rgb(240, 173, 76, maxColorValue = 256)
                    )
                )
            )

            ComplexHeatmap::Heatmap(
                t(datres$m), use_raster = TRUE,
                show_column_names = FALSE, left_annotation = rowAnnot,
                bottom_annotation = colAnnot, name = " "
            )
        })

        output$hintonPlot <- shiny::renderPlot({
            shiny::validate(
                shiny::need(any(!is.na(datres$res$pval)),
                            "No valid p-values")
            )
            tmp <- datres$res
            tmp$batchaff <- factor(tmp$batchaff)
            tmp$condaff <- factor(tmp$condaff)
            j.out <- iSEE::jitterSquarePoints(tmp$condaff, tmp$batchaff)
            summary.data <- j.out$summary
            tmp$jitteredX <- j.out$X
            tmp$jitteredY <- j.out$Y

            ggplot2::ggplot(tmp) +
                ggplot2::geom_tile(
                    ggplot2::aes(x = X, y = Y,
                                 height = 2 * YWidth, width = 2 * XWidth,
                                 group = interaction(X, Y)),
                    summary.data, color = 'black', alpha = 0, size = 0.5
                ) +
                ggplot2::geom_point(
                    aes(x = jitteredX, y = jitteredY, color = padj < 0.05),
                    alpha = 1, tmp, size = 2
                ) +
                ggplot2::labs(x = "Affected by condition",
                              y = "Affected by batch") +
                ggplot2::scale_x_discrete(drop = FALSE) +
                ggplot2::scale_y_discrete(drop = FALSE) +
                ggplot2::scale_color_manual(
                    values = c(`FALSE` = "grey", `TRUE` = "red")
                ) +
                ggplot2::theme_bw() +
                ggplot2::theme(legend.position = 'bottom',
                               legend.text = ggplot2::element_text(size=9),
                               legend.title = ggplot2::element_text(size = 11),
                               legend.box = 'vertical',
                               axis.text = ggplot2::element_text(size = 12),
                               axis.title = ggplot2::element_text(size = 15))

        })

        output$pvalHist <- shiny::renderPlot({
            shiny::validate(
                shiny::need(any(!is.na(datres$res$pval)),
                            "No valid p-values")
            )
            ggplot2::ggplot(
                data = data.frame(pval = datres$res$pval), aes(x = pval)
            ) +
                ggplot2::geom_histogram(bins = 50, fill = "lightgrey") +
                ggplot2::theme_bw() +
                ggplot2::labs(x = "P-value", y = "Count") +
                ggplot2::theme(axis.text = element_text(size = 12),
                               axis.title = element_text(size = 15))
        })

        output$precRecPlot <- shiny::renderPlot({
            shiny::validate(
                shiny::need(any(!is.na(datres$res$pval)),
                            "No valid p-values")
            )
            cbd <- iCOBRA::COBRAData(
                pval = data.frame(mth = datres$res$pval,
                                  row.names = datres$res$feature),
                truth = data.frame(truth = datres$res$condaff,
                                   row.names = datres$res$feature)
            )
            cbd <- iCOBRA::calculate_adjp(cbd)
            cbd <- iCOBRA::calculate_performance(
                cbd, binary_truth = "truth",
                aspects = c("fdrtprcurve", "fdrtpr", "fpr", "tpr")
            )
            cbd <- iCOBRA::prepare_data_for_plot(cbd)
            cowplot::plot_grid(iCOBRA::plot_fpr(cbd),
                               iCOBRA::plot_tpr(cbd))
        })
    }

    # Generate the application
    shiny::shinyApp(ui = ui, server = server)
}
