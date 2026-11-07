#!/usr/bin/env Rscript

# Project-specific wrapper for the current high-coverage dataset.
#
# This file owns the dataset-specific preparation rule. The reusable estimator
# receives only the resulting analysis-ready data object and has no knowledge
# of P_reg.

command_line <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command_line, value = TRUE)
if (!length(script_argument)) {
    stop("Run this file with Rscript.", call. = FALSE)
}
script_path <- normalizePath(
    sub("^--file=", "", script_argument[1L]),
    mustWork = TRUE
)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(repo_root, "scripts", "estimate_pool_effects.R"))

input_file <- file.path(
    repo_root,
    "data",
    "input",
    "birdbath-and-bird-conc0-preg0.0-hicov.txt"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L || (length(args) && args[1L] %in% c("-h", "--help"))) {
    cat(
        paste0(
            "Usage:\n",
            "  Rscript scripts/run_high_coverage_analysis.R ",
            "[BIRDBATH_PREG_THRESHOLD] [OUTPUT_DIRECTORY]\n\n",
            "Defaults: threshold 0.95 and results/final_pool_effect_model.\n"
        )
    )
    quit(status = if (length(args) && args[1L] %in% c("-h", "--help")) 0L else 1L)
}

birdbath_preg_threshold <- if (length(args) >= 1L) {
    suppressWarnings(as.numeric(args[1L]))
} else {
    0.95
}
if (
    length(birdbath_preg_threshold) != 1L ||
        is.na(birdbath_preg_threshold) ||
        !is.finite(birdbath_preg_threshold) ||
        birdbath_preg_threshold < 0 ||
        birdbath_preg_threshold > 1
) {
    stop("BIRDBATH_PREG_THRESHOLD must be a number from 0 through 1.", call. = FALSE)
}

output_dir <- if (length(args) >= 2L) {
    args[2L]
} else {
    file.path(repo_root, "results", "final_pool_effect_model")
}

raw_data <- read.delim(
    input_file,
    header = FALSE,
    sep = "\t",
    na.strings = c("NA", ""),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
if (ncol(raw_data) < 4L) {
    stop(
        "Expected variant, BIRDbath theta, at least one BIRD pool, and final BIRDbath P_reg columns.",
        call. = FALSE
    )
}

birdbath_preg <- suppressWarnings(as.numeric(raw_data[[ncol(raw_data)]]))
invalid_preg <- !is.na(raw_data[[ncol(raw_data)]]) & is.na(birdbath_preg)
invalid_range <- !is.na(birdbath_preg) & (
    !is.finite(birdbath_preg) |
        birdbath_preg < 0 |
        birdbath_preg > 1
)
if (any(invalid_preg | invalid_range)) {
    stop(
        "The final P_reg column must contain numbers from 0 through 1 or NA.",
        call. = FALSE
    )
}

keep <- !is.na(birdbath_preg) & birdbath_preg > birdbath_preg_threshold
ready_data <- raw_data[keep, -ncol(raw_data), drop = FALSE]
if (!nrow(ready_data)) {
    stop(
        "No variants remain after applying BIRDbath P_reg > ",
        birdbath_preg_threshold,
        ".",
        call. = FALSE
    )
}

cat(
    "Prepared", nrow(ready_data), "of", nrow(raw_data),
    "high-coverage variants using BIRDbath P_reg >",
    birdbath_preg_threshold, "\n"
)
estimate_pool_effects(ready_data, output_dir)
cat("Results:", normalizePath(output_dir), "\n")
