#' Confounding explorer
#'
#' @author Charlotte Soneson
#'
#' @param sampleSizes 2x2 numeric matrix giving the number of samples in each
#'   group. Row names must be c('group1', 'group2') and column names must be
#'   c('batch1', 'batch2').
#' @param fracVarCond,fracVarBatch,fracVarUnknown Numeric scalars between 0 and
#'   1. The fraction of variables affected by the condition effect, batch
#'   effect, and 'unknown' effect, respectively.
#' @param condEffectSize,batchEffectSize,unknownEffectSize Numeric scalars.
#'   The condition, batch and 'unknown' effect size, respectively.
#' @param unknownEffectType Character scalar, either 'categorical' or
#'   'continuous', representing the type of 'unknown' effect to add.
#' @param analysisApproach Character scalar. One of 'dontAdjust', 'inclBatch',
#'   'removeBatch', 'removeBatchAccCond'. Determines what model is fit to the
#'   data.
#' @param seed Numeric scalar, the random seed to use when simulating data.
#'
#' @export
#'
#' @importFrom shiny sliderInput radioButtons numericInput fluidRow column
#'   reactiveValues observe validate need renderPlot shinyApp div HTML tags
#'   actionButton plotOutput eventReactive stopApp observeEvent
#' @importFrom ggplot2 ggplot geom_tile aes geom_point labs scale_x_discrete
#'   scale_y_discrete scale_color_manual theme_bw theme element_text
#'   geom_histogram geom_boxplot geom_jitter facet_wrap
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
#' @importFrom dplyr group_by sample_n pull filter left_join select mutate %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom grDevices rgb
#'
#' @return A shinyApp object
#'
#' @examples
#' if (interactive()) {
#'   ConfoundingExplorer()
#' }
#'
ConfoundingExplorer <- function(
    sampleSizes = matrix(rep(5, 4), nrow = 2,
                         dimnames = list(c("group1", "group2"),
                                         c("batch1", "batch2"))),
    fracVarCond = 0.25,
    fracVarBatch = 0.5,
    fracVarUnknown = 0,
    condEffectSize = 3,
    batchEffectSize = 3,
    unknownEffectSize = 0,
    unknownEffectType = "categorical",
    analysisApproach = "dontAdjust",
    seed = 123) {

    .checkInputArguments(sampleSizes = sampleSizes, fracVarCond = fracVarCond,
                         fracVarBatch = fracVarBatch,
                         fracVarUnknown = fracVarUnknown,
                         condEffectSize = condEffectSize,
                         batchEffectSize = batchEffectSize,
                         unknownEffectSize = unknownEffectSize,
                         unknownEffectType = unknownEffectType,
                         analysisApproach = analysisApproach,
                         seed = seed)

    initAppr <-
        switch(analysisApproach,
               dontAdjust = "Don't account for batch effect",
               inclBatch = "Include batch effect in model",
               removeBatch = "Remove batch effect in advance",
               removeBatchAccCond = paste0("Remove batch effect in advance,",
                                           "\naccounting for condition"))

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
                    label = "Group sizes (click outside the table for changes to take effect)",
                    rows = list(
                        names = TRUE,
                        editableNames = FALSE),
                    cols = list(
                        names = TRUE,
                        editableNames = FALSE),
                    value = sampleSizes
                )

            ),
            shinydashboard::menuItem(
                "Set fraction of affected variables",
                startExpanded = TRUE,

                shiny::sliderInput(
                    "fraccond",
                    "Fraction of variables affected by condition:",
                    min = 0, max = 1, value = fracVarCond
                ),
                shiny::sliderInput(
                    "fracbatch",
                    "Fraction of variables affected by batch:",
                    min = 0, max = 1, value = fracVarBatch
                ),
                shiny::sliderInput(
                    "fracunknown",
                    "Fraction of variables affected by unknown variation:",
                    min = 0, max = 1, value = fracVarUnknown
                )
            ),
            shinydashboard::menuItem(
                "Set effect sizes",
                startExpanded = TRUE,

                shiny::sliderInput(
                    "condeffect", "Condition effect size:",
                    min = 0, max = 10, step = 0.1, value = condEffectSize
                ),
                shiny::sliderInput(
                    "batcheffect", "Batch effect size:",
                    min = 0, max = 10, step = 0.1, value = batchEffectSize
                ),
                shiny::sliderInput(
                    "unknowneffect", "Unknown effect size:",
                    min = 0, max = 10, step = 0.1, value = unknownEffectSize
                ),
                shiny::radioButtons(
                    "unknowntype", "Unknown effect type",
                    choices = c("categorical", "continuous"),
                    selected = unknownEffectType, inline = TRUE
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
                selected = initAppr
            ),

            shiny::numericInput(
                "seed", "Random seed", seed, min = 1, max = 1e8
            ),

            shiny::actionButton("datadesc", "Data description"),
            shiny::actionButton("close_app", "Close app")
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
                        shiny::radioButtons("columnorder", "Order features by:",
                                            choices = c("clustering",
                                                        "batch-affected",
                                                        "cond-affected"),
                                            selected = "clustering",
                                            inline = TRUE),
                        shiny::plotOutput("heatmapPlot")
                    )
                )
            ),
            shiny::fluidRow(
                shiny::column(
                    width = 12,
                    shinydashboard::box(
                        id = "example_features",
                        title = shiny::div(
                            "Example features",
                            shiny::HTML("&nbsp;"),
                            shiny::div(id = "help_example_features",
                                       style = "display: inline-block;",
                                       shiny::tags$i(
                                           class = "fa fa-question-circle",
                                           style = "color: lightgrey"
                                       )
                            ),
                            shiny::HTML("&nbsp;"),
                            shiny::HTML("&nbsp;"),
                            shiny::actionButton(inputId = "newExampleFeatures",
                                                label = "New selection")
                        ),
                        width = NULL,
                        shiny::checkboxInput(inputId = "colorExamplesByBatch",
                                             label = "Color samples by batch",
                                             value = TRUE),
                        shiny::plotOutput("exampleFeaturesPlot", height = "23vh")
                    )
                )
            )
        )
    )

    #nocov start
    server <- function(input, output, session) {

        ## Close app
        shiny::observeEvent(input$close_app, {
            shiny::stopApp(returnValue = list(
                res = datres$res, m = datres$m,
                annot = datres$annot
            ))
        })

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
                                 fracunknown = input$fracunknown,
                                 condeffect = input$condeffect,
                                 batcheffect = input$batcheffect,
                                 unknowneffect = input$unknowneffect,
                                 unknowntype = input$unknowntype,
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

            colAnnotDf <- data.frame(batch = datres$res$batchaff,
                                     cond = datres$res$condaff,
                                     row.names = rownames(datres$res))

            if (!any(is.na(datres$res$p.adj))) {
                ## Valid p-values available - add them to heatmap annotation
                colAnnotDf$condsignif <- datres$res$p.adj <= 0.05
            }
            if (any(datres$res$unknownaff)) {
                colAnnotDf$unknown <- datres$res$unknownaff
            }
            colAnnot <- ComplexHeatmap::columnAnnotation(
                df = colAnnotDf,
                col = list(
                    batch = c(`TRUE` = "forestgreen", `FALSE` = "grey85"),
                    cond = c(`TRUE` = "purple", `FALSE` = "grey85"),
                    condsignif = c(`TRUE` = "red", `FALSE` = "grey85"),
                    unknown = c(`TRUE` = "steelblue", `FALSE` = "grey85")
                )
            )

            rowAnnotDf <- data.frame(batch = datres$annot$batch,
                                     cond = datres$annot$cond,
                                     row.names = rownames(datres$annot))
            rowAnnot <- ComplexHeatmap::rowAnnotation(
                df = rowAnnotDf,
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

            ## Determine column order
            if (input$columnorder == "clustering") {
                column_order <- NULL
                cluster_columns <- TRUE
            } else if (input$columnorder == "batch-affected") {
                column_order <- order(paste0(datres$res$batchaff,
                                             datres$res$batchsign),
                                      paste0(datres$res$condaff,
                                             datres$res$condsign),
                                      paste0(datres$res$unknownaff,
                                             datres$res$unknownsign))
                cluster_columns <- FALSE
            } else if (input$columnorder == "cond-affected") {
                column_order <- order(paste0(datres$res$condaff,
                                             datres$res$condsign),
                                      paste0(datres$res$batchaff,
                                             datres$res$batchsign),
                                      paste0(datres$res$unknownaff,
                                             datres$res$unknownsign))
                cluster_columns <- FALSE
            }

            ComplexHeatmap::Heatmap(
                t(datres$m), use_raster = TRUE,
                show_column_names = FALSE, left_annotation = rowAnnot,
                bottom_annotation = colAnnot, name = " ",
                column_order = column_order, cluster_columns = cluster_columns
            )
        })

        ## ----------------------------------------------------------------- ##
        ## Boxplot of example features
        ## ----------------------------------------------------------------- ##
        ## Generate features to plot. Trigger new selection when data is
        ## updated or button is clicked
        feats <- shiny::eventReactive(
            list(input$newExampleFeatures, datres$res), {
                datres$res %>%
                    dplyr::group_by(batchaff, condaff) %>%
                    dplyr::sample_n(1) %>%
                    dplyr::pull(feature)
            }, ignoreInit = FALSE, ignoreNULL = FALSE)

        output$exampleFeaturesPlot <- shiny::renderPlot({
            shiny::validate(
                shiny::need(!is.null(datres$m) & !is.null(datres$res),
                            "Results could not be generated")
            )
            ## Get the data for the selected features
            plotdata <- as.data.frame(datres$m) %>%
                tibble::rownames_to_column("feature") %>%
                dplyr::filter(feature %in% feats()) %>%
                tidyr::gather(key = "sample", value = "value", -feature) %>%
                dplyr::left_join(datres$res %>%
                                     dplyr::select(feature, batchaff,
                                                   condaff), by = "feature") %>%
                dplyr::left_join(datres$annot, by = "sample") %>%
                dplyr::mutate(categ = factor(paste0("batch", batchaff, ".cond",
                                                    condaff),
                                             levels = c("batchFALSE.condFALSE",
                                                        "batchFALSE.condTRUE",
                                                        "batchTRUE.condFALSE",
                                                        "batchTRUE.condTRUE"))) %>%
                dplyr::select(-batchaff, -condaff)
            plt <- ggplot2::ggplot(plotdata, ggplot2::aes(x = cond, y = value)) +
                ggplot2::geom_boxplot(outlier.size = -1) +
                ggplot2::facet_wrap(~ categ, nrow = 1, drop = FALSE) +
                ggplot2::theme_bw()
            if (input$colorExamplesByBatch) {
                plt <- plt +
                    ggplot2::geom_jitter(size = 3, width = 0.2, height = 0,
                                         ggplot2::aes(color = batch)) +
                    ggplot2::scale_color_manual(values = c(
                        B1 = grDevices::rgb(181, 24, 25, maxColorValue = 256),
                        B2 = grDevices::rgb(37, 202, 67, maxColorValue = 256)
                    ))
            } else {
                plt <- plt +
                    ggplot2::geom_jitter(size = 3, width = 0.2, height = 0)
            }
            plt
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
                               legend.text = ggplot2::element_text(size = 9),
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
        ## Help text
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

        shinyjs::onclick("help_example_features", {
            ptour <- data.frame(
                element = "#example_features",
                intro = paste("This plot shows the simulated values of up to ",
                              "four randomly selected features. The panel",
                              "titles indicates whether the feature is",
                              "affected by batch and/or condition effects.",
                              "If possible, one feature from each category is",
                              "randomly sampled. Clicking the",
                              "button samples a new set of features."))
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
    #nocov end

    ## ----------------------------------------------------------------- ##
    ## Generate the app
    ## ----------------------------------------------------------------- ##
    shiny::shinyApp(ui = ui, server = server)
}
