#' @importFrom  methods is
.checkInputArguments <- function(sampleSizes, fracVarCond, fracVarBatch,
                                 fracVarUnknown, condEffectSize,
                                 batchEffectSize, unknownEffectSize,
                                 unknownEffectType, analysisApproach, seed,
                                 addStopButton) {
    if (!methods::is(sampleSizes, "matrix")) {
        stop("sampleSizes must be a matrix")
    }
    if (nrow(sampleSizes) != 2 || ncol(sampleSizes) != 2) {
        stop("sampleSizes must be a 2x2 matrix")
    }
    if (!(all(colnames(sampleSizes) == c("batch1", "batch2")))) {
        stop("The column names of sampleSizes must be batch1, batch2")
    }
    if (!(all(rownames(sampleSizes) == c("group1", "group2")))) {
        stop("The row names of sampleSizes must be group1, group2")
    }
    if (!is.numeric(sampleSizes)) {
        stop("sampleSizes must be numeric")
    }
    if (!all(sampleSizes >= 0)) {
        stop("All elements of sampleSizes must be non-negative")
    }
    sampleSizes <- round(sampleSizes)

    if (!is.numeric(fracVarCond) || length(fracVarCond) > 1 ||
        fracVarCond < 0 || fracVarCond > 1) {
        stop("fracVarCond must be a numeric scalar between 0 and 1")
    }
    if (!is.numeric(fracVarBatch) || length(fracVarBatch) > 1 ||
        fracVarBatch < 0 || fracVarBatch > 1) {
        stop("fracVarBatch must be a numeric scalar between 0 and 1")
    }
    if (!is.numeric(fracVarUnknown) || length(fracVarUnknown) > 1 ||
        fracVarUnknown < 0 || fracVarUnknown > 1) {
        stop("fracVarUnknown must be a numeric scalar between 0 and 1")
    }

    if (!is.numeric(condEffectSize) || length(condEffectSize) > 1 ||
        condEffectSize < 0) {
        stop("condEffectSize must be a non-negative numeric scalar")
    }
    if (!is.numeric(batchEffectSize) || length(batchEffectSize) > 1 ||
        batchEffectSize < 0) {
        stop("batchEffectSize must be a non-negative numeric scalar")
    }
    if (!is.numeric(unknownEffectSize) || length(unknownEffectSize) > 1 ||
        unknownEffectSize < 0) {
        stop("unknownEffectSize must be a non-negative numeric scalar")
    }

    if (!is.character(unknownEffectType) || length(unknownEffectType) != 1 ||
        !(unknownEffectType %in% c("categorical", "continuous"))) {
        stop("unknownEffectType must be one of 'categorical' or 'continuous'")
    }

    if (!is.character(analysisApproach) || length(analysisApproach) != 1 ||
        !(analysisApproach %in% c("dontAdjust", "inclBatch", "removeBatch",
                                  "removeBatchAccCond"))) {
        stop("analysisApproach must be one of 'dontAdjust', 'inclBatch',
             'removeBatch', 'removeBatchAccCond'")
    }

    if (!is.numeric(seed) || length(seed) > 1 || seed <= 0) {
        stop("seed must be a positive numeric scalar")
    }

    if (!is.logical(addStopButton) || length(addStopButton) != 1) {
        stop("addStopButton must be a logical scalar")
    }

}
