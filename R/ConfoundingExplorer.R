#' Confounding explorer
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @importFrom shiny sliderInput radioButtons numericInput fluidRow column
#'   reactiveValues observe validate need renderPlot shinyApp div HTML tags
#' @importFrom ggplot2 ggplot geom_tile aes geom_point labs scale_x_discrete
#'   scale_y_discrete scale_color_manual theme_bw theme element_text
#'   geom_histogram
#' @importFrom iSEE jitterSquarePoints
#' @importFrom iCOBRA COBRAData calculate_adjp calculate_performance
#'   fdrtpr fpr
#' @importFrom cowplot plot_grid
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#'   menuItem dashboardBody box
#' @importFrom ComplexHeatmap columnAnnotation rowAnnotation Heatmap
#' @importFrom shinyMatrix matrixInput
#' @importFrom grDevices rgb
#' @importFrom shinyjs useShinyjs onclick
#' @importFrom rintrojs introjs
#'
#' @return A shinyApp object
#'
#' @examples
#' if (interactive()) {
#'   ConfoundingExplorer()
#' }
#'
ConfoundingExplorer <- function() {
    ui <- shinydashboard::dashboardPage(
        skin = "red",
        shinyjs::useShinyjs(),

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
                                   "\naccounting for condition")),
                selected = "Include batch effect in model"
            ),

            shiny::numericInput(
                "seed", "Random seed", 123, min = 1, max = 1e8
            ),

            shiny::actionButton("datadesc", "Data description")
        ),

        body = shinydashboard::dashboardBody(
            rintrojs::introjsUI(),

            shiny::fluidRow(
                shiny::column(
                    width = 6,
                    shinydashboard::box(
                        id = "statistical_test_results",
                        title = shiny::div(
                            "Statistical test results (Condition effect)",
                            shiny::HTML("&nbsp;"),
                            shiny::div(id = "help_statistical_test_results",
                                       style = "display: inline-block;",
                                       shiny::tags$i(
                                           class = "fa fa-question-circle",
                                           style = "color: lightgrey"
                                       )
                            )
                        ),
                        # title = "Statistical test results (Condition effect)",
                        width = NULL,
                        shiny::plotOutput("hintonPlot")
                    )
                ),
                shiny::column(
                    width = 6,
                    shinydashboard::box(
                        id = "pvalue_histogram",
                        title = shiny::div(
                            "P-value histogram (Condition effect)",
                            shiny::HTML("&nbsp;"),
                            shiny::div(id = "help_pvalue_histogram",
                                       style = "display: inline-block;",
                                       shiny::tags$i(
                                           class = "fa fa-question-circle",
                                           style = "color: lightgrey"
                                       )
                            )
                        ),
                        # title = "P-value histogram (Condition effect)",
                        width = NULL,
                        shiny::plotOutput("pvalHist")
                    )
                )
            ),
            shiny::fluidRow(
                shiny::column(
                    width = 12,
                    shinydashboard::box(
                        id = "summary",
                        title = shiny::div(
                            "Summary",
                            shiny::HTML("&nbsp;"),
                            shiny::div(id = "help_summary",
                                       style = "display: inline-block;",
                                       shiny::tags$i(
                                           class = "fa fa-question-circle",
                                           style = "color: lightgrey"
                                       )
                            )
                        ),
                        # title = "Summary",
                        width = NULL,
                        shiny::textOutput("summaryRates")
                    )
                )
            ),
            shiny::fluidRow(
                shiny::column(
                    width = 12,
                    shinydashboard::box(
                        id = "data_heatmap",
                        title = shiny::div(
                            "Data heatmap",
                            shiny::HTML("&nbsp;"),
                            shiny::div(id = "help_data_heatmap",
                                       style = "display: inline-block;",
                                       shiny::tags$i(
                                           class = "fa fa-question-circle",
                                           style = "color: lightgrey"
                                       )
                            )
                        ),
                        # title = "Data heatmap",
                        width = NULL,
                        shiny::plotOutput("heatmapPlot")
                    )
                )
            )
        )
    )

    server <- function(input, output, session) {

        ## Display data generation description
        shiny::observeEvent(input$datadesc, {
            shiny::showModal(shiny::modalDialog(
                title = "Data generation",
                shiny::includeHTML(system.file(
                    "extdata/dataDescription.html",
                    package = "ConfoundingExplorer")
                ),
                # shiny::renderUI(shiny::HTML(readLines(
                #     system.file("extdata/dataDescription.html",
                #                 package = "ConfoundingExplorer")))),
                # shiny::HTML(readLines(
                #     system.file("extdata/dataDescription.html",
                #                 package = "ConfoundingExplorer"))),
                # shiny::includeMarkdown(),
                easyClose = TRUE
            ))
        }, ignoreNULL = TRUE, ignoreInit = TRUE)

        ## ----------------------------------------------------------------- ##
        ## Data generation
        ## ----------------------------------------------------------------- ##
        datres <- shiny::reactiveValues()
        shiny::observe({
            shiny::validate(
                shiny::need(is.numeric(input$seed),
                            "Seed must be a number")
            )
            tmp <- .generateData(nb1c1 = input$sampleSizeTable[1, 1],
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

        ## ----------------------------------------------------------------- ##
        ## Data heatmap
        ## ----------------------------------------------------------------- ##
        output$heatmapPlot <- shiny::renderPlot({
            shiny::validate(
                shiny::need(!is.null(datres$m) & !is.null(datres$res),
                            "Results could not be generated")
            )
            colAnnot <- ComplexHeatmap::columnAnnotation(
                batch = datres$res$batchaff,
                cond = datres$res$condaff,
                col = list(
                    batch = c(`TRUE` = "forestgreen", `FALSE` = "grey85"),
                    cond = c(`TRUE` = "purple", `FALSE` = "grey85")
                )
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

        ## ----------------------------------------------------------------- ##
        ## Hinton plot of variables (batch/condition affected, adjp<0.05)
        ## ----------------------------------------------------------------- ##
        output$hintonPlot <- shiny::renderPlot({
            shiny::validate(
                shiny::need(any(!is.na(datres$res$p.val)),
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
                    aes(x = jitteredX, y = jitteredY, color = p.adj < 0.05),
                    alpha = 1, tmp, size = 2
                ) +
                ggplot2::labs(x = "Variable affected by condition",
                              y = "Variable affected by batch") +
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

        ## ----------------------------------------------------------------- ##
        ## P-value histogram
        ## ----------------------------------------------------------------- ##
        output$pvalHist <- shiny::renderPlot({
            shiny::validate(
                shiny::need(any(!is.na(datres$res$p.val)),
                            "No valid p-values")
            )
            ggplot2::ggplot(
                data = data.frame(p.val = datres$res$p.val), aes(x = p.val)
            ) +
                ggplot2::geom_histogram(bins = 50, fill = "lightgrey") +
                ggplot2::theme_bw() +
                ggplot2::labs(x = "P-value", y = "Count") +
                ggplot2::theme(axis.text = element_text(size = 12),
                               axis.title = element_text(size = 15))
        })

        ## ----------------------------------------------------------------- ##
        ## Summary (FDP, TPR, FPR)
        ## ----------------------------------------------------------------- ##
        output$summaryRates <- shiny::renderText({
            shiny::validate(
                shiny::need(any(!is.na(datres$res$p.val)),
                            "No valid p-values")
            )
            cbd <- iCOBRA::COBRAData(
                pval = data.frame(mth = datres$res$p.val,
                                  row.names = datres$res$feature),
                padj = data.frame(mth = datres$res$p.adj,
                                  row.names = datres$res$feature),
                truth = data.frame(truth = datres$res$condaff,
                                   row.names = datres$res$feature)
            )
            suppressMessages({
                cbd <- iCOBRA::calculate_performance(
                    cbd, binary_truth = "truth", thrs = 0.05,
                    aspects = c("fpr", "fdrtpr")
                )
            })
            paste0("At adj.p threshold = 0.05: TPR = ",
                   signif(iCOBRA::fdrtpr(cbd)$TPR, digits = 3), ", FDP = ",
                   signif(iCOBRA::fdrtpr(cbd)$FDR, digits = 3), ", FPR = ",
                   signif(iCOBRA::fpr(cbd)$FPR, digits = 3))
        })

        ## ----------------------------------------------------------------- ##
        ## Help textx
        ## ----------------------------------------------------------------- ##
        shinyjs::onclick("help_statistical_test_results", {
            ptour <- data.frame(
                element = "#statistical_test_results",
                intro = paste("This plot shows the number of variables in",
                              "each category (affected by batch and/or",
                              "condition). The sizes of the rectangles are",
                              "proportional to the number of variables in",
                              "the respective categories. Each dot represents",
                              "a variable. Variables with an adjusted p-value",
                              "below 0.05 (testing the null hypothesis that",
                              "there is no difference between conditions)",
                              "are colored in red."))
            rintrojs::introjs(session, options = list(steps = ptour))
        })

        shinyjs::onclick("help_pvalue_histogram", {
            ptour <- data.frame(
                element = "#pvalue_histogram",
                intro = paste("This plot shows a histogram of the nominal",
                              "pvalues (testing the null hypothesis that",
                              "there is no difference between conditions)."))
            rintrojs::introjs(session, options = list(steps = ptour))
        })

        shinyjs::onclick("help_summary", {
            ptour <- data.frame(
                element = "#summary",
                intro = paste("This panel summarizes the true positive rate",
                              "(TPR, power, the fraction of the truly",
                              "differential variables that are detected as",
                              "such), the false discovery proportion (FDP,",
                              "the fraction of the significant variables that",
                              "are false positives) and the false positive",
                              "rate (FPR, the fraction of the truly",
                              "non-differential variables that are called",
                              "significant) at an imposed adjusted p-value",
                              "threshold of 0.05."))
            rintrojs::introjs(session, options = list(steps = ptour))
        })

        shinyjs::onclick("help_data_heatmap", {
            ptour <- data.frame(
                element = "#data_heatmap",
                intro = paste("This plot provides an overview of the simulated",
                              "data set. Each row is a sample (the",
                              "associated condition and batch are indicated),",
                              "and each column is a variable (the colored bars",
                              "underneath the plot indicates which variables",
                              "are truly affected by the condition and batch,",
                              "respectively)."))
            rintrojs::introjs(session, options = list(steps = ptour))
        })

    }

    ## ----------------------------------------------------------------- ##
    ## Generate the app
    ## ----------------------------------------------------------------- ##
    shiny::shinyApp(ui = ui, server = server)
}
