parse_optional_preg_threshold <- function(value) {
    if (is.null(value) || !length(value) || identical(value, "")) {
        return(NULL)
    }

    threshold <- suppressWarnings(as.numeric(value))
    if (
        length(threshold) != 1L ||
            is.na(threshold) ||
            !is.finite(threshold) ||
            threshold < 0 ||
            threshold > 1
    ) {
        stop(
            "BIRDBATH_PREG_THRESHOLD must be a number from 0 through 1.",
            call. = FALSE
        )
    }
    threshold
}


read_birdbath_wide_input <- function(
    input_file,
    birdbath_preg_threshold = NULL
) {
    if (!file.exists(input_file)) {
        stop("Input file does not exist: ", input_file, call. = FALSE)
    }

    wide <- read.delim(
        input_file,
        header = FALSE,
        sep = "\t",
        na.strings = c("NA", ""),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    variants_before_filter <- nrow(wide)
    preg_available <- !is.null(birdbath_preg_threshold)

    if (preg_available) {
        if (ncol(wide) < 4L) {
            stop(
                "A P_reg-aware input requires variant, BIRDbath theta, at least one BIRD pool, and final BIRDbath P_reg columns.",
                call. = FALSE
            )
        }

        original <- wide[[ncol(wide)]]
        birdbath_preg <- suppressWarnings(as.numeric(original))
        invalid_text <- !is.na(original) & is.na(birdbath_preg)
        invalid_range <- !is.na(birdbath_preg) & (
            !is.finite(birdbath_preg) |
                birdbath_preg < 0 |
                birdbath_preg > 1
        )
        if (any(invalid_text | invalid_range)) {
            stop(
                "The final BIRDbath P_reg column must contain numbers from 0 through 1 or NA.",
                call. = FALSE
            )
        }

        retained <- !is.na(birdbath_preg) &
            birdbath_preg > birdbath_preg_threshold
        wide <- wide[retained, -ncol(wide), drop = FALSE]
        if (!nrow(wide)) {
            stop(
                "No variants remain after applying BIRDbath P_reg > ",
                birdbath_preg_threshold,
                ".",
                call. = FALSE
            )
        }
    }

    list(
        wide = wide,
        birdbath_preg_available = preg_available,
        birdbath_preg_threshold = birdbath_preg_threshold,
        variants_before_preg_filter = variants_before_filter,
        variants_after_preg_filter = nrow(wide),
        variants_excluded_by_preg = variants_before_filter - nrow(wide)
    )
}


preg_filter_summary_lines <- function(data) {
    if (isTRUE(data$birdbath_preg_available)) {
        c(
            sprintf(
                "BIRDbath P_reg filter: P_reg > %s (strict)",
                format(data$birdbath_preg_threshold, trim = TRUE)
            ),
            paste(
                "Variants before BIRDbath P_reg filtering:",
                data$variants_before_preg_filter
            ),
            paste(
                "Variants after BIRDbath P_reg filtering:",
                data$variants_after_preg_filter
            ),
            paste(
                "Variants excluded by BIRDbath P_reg filtering:",
                data$variants_excluded_by_preg
            )
        )
    } else {
        "BIRDbath P_reg filter: not applied; all columns after BIRDbath theta were interpreted as BIRD pools"
    }
}
