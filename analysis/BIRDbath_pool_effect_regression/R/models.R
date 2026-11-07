cluster_robust_vcov <- function(fit, cluster) {
    design <- model.matrix(fit)
    residual <- residuals(fit)

    if (length(cluster) != nrow(design)) {
        stop("Cluster vector and model matrix have different lengths.")
    }

    cluster <- factor(cluster)
    n_clusters <- nlevels(cluster)
    n_observations <- nrow(design)
    n_parameters <- fit$rank

    if (n_clusters <= 1L) {
        stop("At least two variant clusters are required.")
    }

    cluster_scores <- rowsum(
        design * residual,
        group = cluster,
        reorder = FALSE
    )

    bread <- solve(crossprod(design))
    meat <- crossprod(cluster_scores)
    correction <-
        (n_clusters / (n_clusters - 1)) *
        ((n_observations - 1) / (n_observations - n_parameters))

    bread %*% meat %*% bread * correction
}


fit_pool_effect_models <- function(model_data) {
    list(
        fixed_gamma_1 = lm(
            log2_theta_bird ~
                0 + pool + offset(log2_theta_birdbath),
            data = model_data
        ),
        free_gamma = lm(
            log2_theta_bird ~
                0 + pool + log2_theta_birdbath,
            data = model_data
        )
    )
}


calculate_fit_statistics <- function(fit, observed, model_name) {
    residual <- residuals(fit)
    rss <- sum(residual^2)
    tss <- sum((observed - mean(observed))^2)

    data.frame(
        model = model_name,
        n = nobs(fit),
        parameters = fit$rank,
        residual_df = df.residual(fit),
        residual_sd = sqrt(rss / df.residual(fit)),
        rmse = sqrt(mean(residual^2)),
        mae = mean(abs(residual)),
        r_squared = 1 - rss / tss,
        adjusted_r_squared = 1 -
            (rss / df.residual(fit)) /
            (tss / (length(observed) - 1)),
        AIC = AIC(fit),
        BIC = BIC(fit),
        stringsAsFactors = FALSE
    )
}


extract_pool_effects <- function(
    fit,
    cluster_vcov,
    model_name,
    pool_counts,
    cluster,
    confidence_level = 0.95,
    minimum_pool_n_for_warning = 30L
) {
    coefficient_matrix <- summary(fit)$coefficients
    coefficient_names <- rownames(coefficient_matrix)
    pool_terms <- coefficient_names[grepl("^pool", coefficient_names)]

    estimate <- coefficient_matrix[pool_terms, "Estimate"]
    ols_se <- coefficient_matrix[pool_terms, "Std. Error"]
    ols_p_value <- coefficient_matrix[pool_terms, "Pr(>|t|)"]
    robust_se <- sqrt(diag(cluster_vcov))[pool_terms]
    robust_df <- length(unique(cluster)) - 1L
    robust_t <- estimate / robust_se
    robust_p_value <- 2 * pt(
        abs(robust_t),
        df = robust_df,
        lower.tail = FALSE
    )
    critical_value <- qt(
        1 - (1 - confidence_level) / 2,
        df = robust_df
    )

    result <- data.frame(
        pool = sub("^pool", "", pool_terms),
        model = model_name,
        beta = unname(estimate),
        cluster_robust_se = unname(robust_se),
        t_value = unname(robust_t),
        p_value = unname(robust_p_value),
        ols_se = unname(ols_se),
        ols_p_value = unname(ols_p_value),
        ci_lower = unname(estimate - critical_value * robust_se),
        ci_upper = unname(estimate + critical_value * robust_se),
        multiplier_at_birdbath_theta_1 = 2^unname(estimate),
        multiplier_ci_lower =
            2^unname(estimate - critical_value * robust_se),
        multiplier_ci_upper =
            2^unname(estimate + critical_value * robust_se),
        stringsAsFactors = FALSE
    )

    result$p_value_fdr <- p.adjust(result$p_value, method = "BH")
    count_index <- match(result$pool, pool_counts$pool)
    result <- cbind(
        result,
        pool_counts[count_index, setdiff(names(pool_counts), "pool"), drop = FALSE]
    )
    result$direction <- ifelse(
        result$beta > 0,
        "positive",
        ifelse(result$beta < 0, "negative", "zero")
    )
    result$warning_status <- ifelse(
        result$n < minimum_pool_n_for_warning,
        "insufficient_data",
        ifelse(
            result$p_value_fdr < 0.05,
            "pool_effect_detected",
            ifelse(
                result$p_value < 0.05,
                "nominal_only",
                "no_pool_effect_detected"
            )
        )
    )

    result
}


assemble_pool_effect_matrix <- function(
    pool_counts,
    fixed_effects,
    free_effects
) {
    common <- pool_counts[c("pool", "n", "n_variants", "data_status")]
    fields <- c(
        "beta",
        "cluster_robust_se",
        "p_value",
        "p_value_fdr",
        "ci_lower",
        "ci_upper",
        "multiplier_at_birdbath_theta_1",
        "multiplier_ci_lower",
        "multiplier_ci_upper",
        "warning_status"
    )

    fixed <- fixed_effects[c("pool", fields)]
    free <- free_effects[c("pool", fields)]
    names(fixed)[-1L] <- paste0("fixed_", names(fixed)[-1L])
    names(free)[-1L] <- paste0("free_", names(free)[-1L])

    result <- merge(common, fixed, by = "pool", all.x = TRUE, sort = FALSE)
    result <- merge(result, free, by = "pool", all.x = TRUE, sort = FALSE)
    result[match(pool_counts$pool, result$pool), , drop = FALSE]
}


summarize_model_fits <- function(
    fits,
    vcov_free,
    model_data,
    confidence_level = 0.95
) {
    observed <- model_data$log2_theta_bird
    fit_statistics <- rbind(
        calculate_fit_statistics(
            fits$fixed_gamma_1,
            observed,
            "fixed_gamma_1"
        ),
        calculate_fit_statistics(
            fits$free_gamma,
            observed,
            "free_gamma"
        )
    )

    gamma_hat <- unname(coef(fits$free_gamma)["log2_theta_birdbath"])
    gamma_se <- sqrt(
        vcov_free["log2_theta_birdbath", "log2_theta_birdbath"]
    )
    gamma_df <- length(unique(model_data$variant)) - 1L
    critical <- qt(1 - (1 - confidence_level) / 2, df = gamma_df)
    gamma_t <- (gamma_hat - 1) / gamma_se
    gamma_p <- 2 * pt(abs(gamma_t), df = gamma_df, lower.tail = FALSE)

    fixed_rss <- sum(residuals(fits$fixed_gamma_1)^2)
    free_rss <- sum(residuals(fits$free_gamma)^2)
    total_ss <- sum((observed - mean(observed))^2)
    df_added <- fits$free_gamma$rank - fits$fixed_gamma_1$rank
    nested_f <- ((fixed_rss - free_rss) / df_added) /
        (free_rss / df.residual(fits$free_gamma))

    list(
        fit_statistics = fit_statistics,
        gamma = data.frame(
            gamma = gamma_hat,
            cluster_robust_se = gamma_se,
            ci_lower = gamma_hat - critical * gamma_se,
            ci_upper = gamma_hat + critical * gamma_se,
            p_for_gamma_equals_1 = gamma_p,
            twofold_bird_change = 2^gamma_hat,
            stringsAsFactors = FALSE
        ),
        comparison = data.frame(
            parameters_added = df_added,
            delta_r_squared = (fixed_rss - free_rss) / total_ss,
            partial_r_squared = (fixed_rss - free_rss) / fixed_rss,
            rmse_reduction_percent = 100 * (
                1 -
                    fit_statistics$rmse[fit_statistics$model == "free_gamma"] /
                    fit_statistics$rmse[fit_statistics$model == "fixed_gamma_1"]
            ),
            F = nested_f,
            p_value = pf(
                nested_f,
                df1 = df_added,
                df2 = df.residual(fits$free_gamma),
                lower.tail = FALSE
            ),
            stringsAsFactors = FALSE
        )
    )
}
