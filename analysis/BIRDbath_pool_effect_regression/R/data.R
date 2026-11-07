read_numeric_column <- function(x, column_name) {
    converted <- suppressWarnings(as.numeric(x))
    invalid <- !is.na(x) & is.na(converted)

    if (any(invalid)) {
        stop(
            sprintf(
                "Column '%s' contains non-numeric values. The input must be a headerless tab-delimited file.",
                column_name
            ),
            call. = FALSE
        )
    }

    converted
}


read_pool_effect_data <- function(
    input_file,
    minimum_pool_n_for_warning = 30L,
    birdbath_preg_threshold = NULL
) {
    input_is_object <- is.data.frame(input_file) || is.matrix(input_file)

    if (input_is_object) {
        if (!is.null(birdbath_preg_threshold)) {
            stop(
                "P_reg filtering is only supported for legacy file-based workflows. Supply a ready-to-analyze data object instead.",
                call. = FALSE
            )
        }
        input <- list(
            wide = as.data.frame(
                input_file,
                stringsAsFactors = FALSE,
                check.names = FALSE
            ),
            birdbath_preg_available = FALSE,
            birdbath_preg_threshold = NULL,
            variants_before_preg_filter = nrow(input_file),
            variants_after_preg_filter = nrow(input_file),
            variants_excluded_by_preg = 0L
        )
        input_label <- "<in-memory data>"
    } else {
        if (!is.character(input_file) || length(input_file) != 1L) {
            stop(
                "Input must be a ready-to-analyze file path, data frame, or matrix.",
                call. = FALSE
            )
        }
        input <- read_birdbath_wide_input(
            input_file,
            birdbath_preg_threshold
        )
        input_label <- normalizePath(input_file)
    }
    wide <- input$wide

    if (ncol(wide) < 3L) {
        stop(
            "Expected at least three columns: variant, BIRDbath theta, and one or more BIRD pool columns.",
            call. = FALSE
        )
    }

    pool_levels <- paste0("pool", seq_len(ncol(wide) - 2L))
    colnames(wide) <- c("variant", "theta_birdbath", pool_levels)
    wide$variant <- as.character(wide$variant)
    wide$theta_birdbath <- read_numeric_column(
        wide$theta_birdbath,
        "theta_birdbath"
    )

    for (pool_name in pool_levels) {
        wide[[pool_name]] <- read_numeric_column(
            wide[[pool_name]],
            pool_name
        )
    }

    if (anyDuplicated(wide$variant)) {
        warning(
            "The input contains duplicated variant identifiers. Cross-validation will keep identical identifiers in the same fold.",
            call. = FALSE
        )
    }

    n_variants <- nrow(wide)
    n_pools <- length(pool_levels)

    raw_long <- data.frame(
        variant = rep(wide$variant, times = n_pools),
        theta_birdbath = rep(wide$theta_birdbath, times = n_pools),
        pool = factor(
            rep(pool_levels, each = n_variants),
            levels = pool_levels
        ),
        theta_bird = unlist(
            wide[pool_levels],
            use.names = FALSE
        ),
        stringsAsFactors = FALSE
    )

    invalid_value <- (
        !is.na(raw_long$theta_birdbath) &
            (!is.finite(raw_long$theta_birdbath) |
                raw_long$theta_birdbath <= 0)
    ) | (
        !is.na(raw_long$theta_bird) &
            (!is.finite(raw_long$theta_bird) |
                raw_long$theta_bird <= 0)
    )

    usable <- !is.na(raw_long$theta_birdbath) &
        !is.na(raw_long$theta_bird) &
        is.finite(raw_long$theta_birdbath) &
        is.finite(raw_long$theta_bird) &
        raw_long$theta_birdbath > 0 &
        raw_long$theta_bird > 0

    model_data <- raw_long[usable, , drop = FALSE]
    model_data$pool <- factor(
        model_data$pool,
        levels = pool_levels
    )
    model_data$log2_theta_birdbath <- log2(model_data$theta_birdbath)
    model_data$log2_theta_bird <- log2(model_data$theta_bird)
    model_data$log2_ratio <-
        model_data$log2_theta_bird - model_data$log2_theta_birdbath

    if (nrow(model_data) == 0L) {
        stop("No usable positive theta pairs remain after validation.", call. = FALSE)
    }

    pool_missingness <- do.call(
        rbind,
        lapply(pool_levels, function(pool_name) {
            values <- raw_long$theta_bird[raw_long$pool == pool_name]
            n_observed <- sum(!is.na(values))

            data.frame(
                pool = pool_name,
                n_total = length(values),
                n_observed = n_observed,
                n_missing = sum(is.na(values)),
                proportion_observed = n_observed / length(values),
                stringsAsFactors = FALSE
            )
        })
    )

    pool_counts <- do.call(
        rbind,
        lapply(pool_levels, function(pool_name) {
            this_pool <- model_data[model_data$pool == pool_name, , drop = FALSE]
            n <- nrow(this_pool)

            data.frame(
                pool = pool_name,
                n = n,
                n_variants = length(unique(this_pool$variant)),
                n_negative_birdbath = sum(this_pool$log2_theta_birdbath < 0),
                n_positive_birdbath = sum(this_pool$log2_theta_birdbath > 0),
                median_log2_birdbath = if (n) {
                    median(this_pool$log2_theta_birdbath)
                } else {
                    NA_real_
                },
                median_log2_bird = if (n) {
                    median(this_pool$log2_theta_bird)
                } else {
                    NA_real_
                },
                data_status = if (n < minimum_pool_n_for_warning) {
                    "insufficient_data"
                } else {
                    "sufficient_for_screening"
                },
                stringsAsFactors = FALSE
            )
        })
    )

    empty_pools <- pool_counts$pool[pool_counts$n == 0L]
    if (length(empty_pools)) {
        stop(
            "No usable observations were found for: ",
            paste(empty_pools, collapse = ", "),
            ". Remove empty columns or provide data for every pool.",
            call. = FALSE
        )
    }

    list(
        input_file = input_label,
        wide = wide,
        model_data = model_data,
        pool_levels = pool_levels,
        pool_counts = pool_counts,
        pool_missingness = pool_missingness,
        invalid_value_count = sum(invalid_value),
        wide_variant_count = nrow(wide),
        birdbath_preg_available = input$birdbath_preg_available,
        birdbath_preg_threshold = input$birdbath_preg_threshold,
        variants_before_preg_filter = input$variants_before_preg_filter,
        variants_after_preg_filter = input$variants_after_preg_filter,
        variants_excluded_by_preg = input$variants_excluded_by_preg
    )
}
