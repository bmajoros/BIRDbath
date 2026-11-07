#!/usr/bin/env Rscript

usage <- function() {
    cat(
        paste0(
            "Usage:\n",
            "  Rscript scripts/estimate_pool_effects.R ",
            "INPUT_FILE OUTPUT_DIRECTORY\n\n",
            "INPUT_FILE must already be ready for analysis: variant ID, ",
            "BIRDbath theta, then one BIRD-theta column per pool.\n",
            "This estimator does not filter rows or interpret a P_reg column.\n"
        )
    )
}

source_path <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
was_sourced <- !is.null(source_path)
if (was_sourced) {
    script_path <- normalizePath(source_path, mustWork = TRUE)
} else {
    command_line <- commandArgs(trailingOnly = FALSE)
    script_argument <- grep("^--file=", command_line, value = TRUE)
    if (!length(script_argument)) {
        stop("Run this file with Rscript or source it from R.", call. = FALSE)
    }
    script_path <- normalizePath(
        sub("^--file=", "", script_argument[1L]),
        mustWork = TRUE
    )
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "input.R"))
source(file.path(repo_root, "R", "data.R"))
source(file.path(repo_root, "R", "models.R"))
source(file.path(repo_root, "R", "final_pool_effect_model.R"))

run_estimate_pool_effects_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
    if (length(args) != 2L || args[1L] %in% c("-h", "--help")) {
        usage()
        quit(
            status = if (
                length(args) && args[1L] %in% c("-h", "--help")
            ) 0L else 1L
        )
    }

    input_file <- args[1L]
    output_dir <- args[2L]

    message("Reading analysis-ready input and fitting the pool-effect model...")
    result <- estimate_pool_effects(input_file, output_dir)
    analysis <- result$analysis
    input_data <- result$input_data

    cat("\nPool-effect endpoint complete.\n")
    cat("Usable observations:", nrow(input_data$model_data), "\n")
    cat("Alpha:", sprintf("%.6f", analysis$model_summary$alpha), "\n")
    cat("Gamma:", sprintf("%.6f", analysis$model_summary$gamma), "\n")
    cat("Results:", normalizePath(output_dir), "\n")
}

if (!was_sourced) {
    run_estimate_pool_effects_cli()
}
