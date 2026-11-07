fit_final_pool_effect_model <- function(
    model_data,
    pool_levels,
    confidence_level = 0.95
) {
    fit_data <- data.frame(
        variant = model_data$variant,
        pool = factor(model_data$pool, levels = pool_levels),
        x = model_data$log2_theta_birdbath,
        y = model_data$log2_theta_bird
    )

    # Fitting one intercept per pool is equivalent to
    #
    #   y = alpha + gamma*x + beta_pool + error,
    #   sum(beta_pool) = 0.
    #
    # For reporting, alpha is the mean of the fitted pool intercepts and each
    # beta_pool is that pool's deviation from the mean. Retaining alpha avoids
    # forcing the average-pool relationship through x = y = 0.
    fit <- lm(y ~ 0 + pool + x, data = fit_data)
    robust_vcov <- cluster_robust_vcov(fit, fit_data$variant)

    coefficient <- coef(fit)
    pool_terms <- names(coefficient)[grepl("^pool", names(coefficient))]
    pool_intercepts <- coefficient[pool_terms]
    n_pools <- length(pool_intercepts)
    center <- diag(n_pools) - matrix(1 / n_pools, n_pools, n_pools)
    beta <- as.vector(center %*% pool_intercepts)
    names(beta) <- pool_levels

    pool_vcov <- robust_vcov[pool_terms, pool_terms, drop = FALSE]
    beta_vcov <- center %*% pool_vcov %*% t(center)
    beta_se <- sqrt(diag(beta_vcov))
    cluster_df <- length(unique(fit_data$variant)) - 1L
    critical <- qt(1 - (1 - confidence_level) / 2, df = cluster_df)
    beta_t <- beta / beta_se
    beta_p <- 2 * pt(abs(beta_t), df = cluster_df, lower.tail = FALSE)

    counts <- do.call(
        rbind,
        lapply(pool_levels, function(pool_name) {
            rows <- fit_data$pool == pool_name
            data.frame(
                pool = pool_name,
                observations = sum(rows),
                variants = length(unique(fit_data$variant[rows])),
                stringsAsFactors = FALSE
            )
        })
    )

    pool_effects <- data.frame(
        pool = pool_levels,
        beta_log2 = beta,
        cluster_robust_se = beta_se,
        ci_lower = beta - critical * beta_se,
        ci_upper = beta + critical * beta_se,
        p_value = beta_p,
        p_value_fdr = p.adjust(beta_p, method = "BH"),
        theta_multiplier_vs_average_pool = 2^beta,
        multiplier_ci_lower = 2^(beta - critical * beta_se),
        multiplier_ci_upper = 2^(beta + critical * beta_se),
        observations = counts$observations,
        variants = counts$variants,
        stringsAsFactors = FALSE
    )

    alpha_weights <- rep(1 / n_pools, n_pools)
    alpha <- mean(pool_intercepts)
    alpha_se <- sqrt(as.numeric(
        t(alpha_weights) %*% pool_vcov %*% alpha_weights
    ))
    gamma <- unname(coefficient["x"])
    gamma_se <- sqrt(robust_vcov["x", "x"])
    residual <- residuals(fit)
    total_ss <- sum((fit_data$y - mean(fit_data$y))^2)

    model_summary <- data.frame(
        observations = nrow(fit_data),
        variants = length(unique(fit_data$variant)),
        pools = n_pools,
        alpha = alpha,
        alpha_cluster_robust_se = alpha_se,
        alpha_ci_lower = alpha - critical * alpha_se,
        alpha_ci_upper = alpha + critical * alpha_se,
        gamma = gamma,
        gamma_cluster_robust_se = gamma_se,
        gamma_ci_lower = gamma - critical * gamma_se,
        gamma_ci_upper = gamma + critical * gamma_se,
        residual_sd = sqrt(sum(residual^2) / df.residual(fit)),
        r_squared = 1 - sum(residual^2) / total_ss,
        stringsAsFactors = FALSE
    )

    list(
        fit = fit,
        fit_data = fit_data,
        robust_vcov = robust_vcov,
        pool_intercepts = setNames(as.numeric(pool_intercepts), pool_levels),
        pool_effects = pool_effects,
        model_summary = model_summary
    )
}


estimate_pool_effects <- function(
    input,
    output_dir,
    confidence_level = 0.95
) {
    # `input` is deliberately assumed to be analysis-ready. It may be either
    # a headerless wide-data file or an in-memory data frame/matrix with
    # variant, BIRDbath theta, and one column per BIRD pool. Dataset-specific
    # filtering belongs in an upstream wrapper, not in this estimator.
    input_data <- read_pool_effect_data(input)
    analysis <- fit_final_pool_effect_model(
        input_data$model_data,
        input_data$pool_levels,
        confidence_level = confidence_level
    )
    write_final_pool_effect_outputs(analysis, input_data, output_dir)

    invisible(list(
        analysis = analysis,
        input_data = input_data,
        output_dir = normalizePath(output_dir)
    ))
}


make_final_pool_effect_plot <- function(pool_effects) {
    plot_data <- pool_effects
    plot_data$pool <- factor(
        plot_data$pool,
        levels = rev(plot_data$pool)
    )

    ggplot2::ggplot(plot_data, ggplot2::aes(y = pool, x = beta_log2)) +
        ggplot2::geom_vline(
            xintercept = 0,
            linetype = 2,
            colour = "grey45"
        ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
            width = 0.25,
            orientation = "y"
        ) +
        ggplot2::geom_point(
            shape = 21,
            size = 3,
            fill = "#356A8A"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::labs(
            x = "Pool effect relative to the average pool, log2 theta units",
            y = "Pool",
            title = "Estimated pool-associated effects",
            subtitle = "Points are beta coefficients; error bars are variant-clustered 95% confidence intervals"
        )
}


make_final_scatter_plot <- function(analysis, pool_levels) {
    fitted_lines <- do.call(
        rbind,
        lapply(pool_levels, function(pool_name) {
            this_pool <- analysis$fit_data[
                analysis$fit_data$pool == pool_name,
                ,
                drop = FALSE
            ]
            x_grid <- seq(
                min(this_pool$x),
                max(this_pool$x),
                length.out = 100L
            )
            data.frame(
                pool = factor(pool_name, levels = pool_levels),
                x = x_grid,
                predicted = analysis$pool_intercepts[pool_name] +
                    analysis$model_summary$gamma * x_grid
            )
        })
    )

    ggplot2::ggplot(
        analysis$fit_data,
        ggplot2::aes(x = x, y = y)
    ) +
        ggplot2::geom_point(
            alpha = 0.20,
            size = 0.55,
            colour = "grey30"
        ) +
        ggplot2::geom_abline(
            intercept = 0,
            slope = 1,
            linetype = 3,
            colour = "grey50"
        ) +
        ggplot2::geom_line(
            data = fitted_lines,
            ggplot2::aes(x = x, y = predicted),
            colour = "#C62828",
            linewidth = 0.8
        ) +
        ggplot2::facet_wrap(~pool, ncol = 5) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::labs(
            x = "log2 BIRDbath theta",
            y = "log2 BIRD theta",
            title = "BIRD versus BIRDbath by pool",
            subtitle = "Dotted line is y = x; red line is the fitted common-slope model"
        )
}


write_final_model_summary <- function(
    input_data,
    model_summary,
    output_path
) {
    lines <- c(
        "FINAL POOL-EFFECT MODEL",
        "",
        paste("Input:", input_data$input_file),
        paste("Usable observations:", model_summary$observations),
        paste("Unique variants:", model_summary$variants),
        paste("Pools:", model_summary$pools),
        "",
        "MODEL",
        "log2(theta_BIRD[v,p]) = alpha + gamma*log2(theta_BIRDbath[v]) + beta_pool[p] + error[v,p]",
        "Pool effects beta_pool are constrained to sum to zero.",
        "Alpha is the fitted average-pool intercept; retaining it avoids forcing the average-pool relationship through the origin.",
        "Gamma is the common BIRDbath slope, and beta_pool is each pool's deviation from the average pool.",
        "Confidence intervals use variant-clustered standard errors.",
        "",
        "ESTIMATES",
        sprintf(
            "Alpha: %.6f (95%% CI %.6f to %.6f)",
            model_summary$alpha,
            model_summary$alpha_ci_lower,
            model_summary$alpha_ci_upper
        ),
        sprintf(
            "Gamma: %.6f (95%% CI %.6f to %.6f)",
            model_summary$gamma,
            model_summary$gamma_ci_lower,
            model_summary$gamma_ci_upper
        ),
        sprintf("Residual SD: %.6f", model_summary$residual_sd),
        sprintf("R-squared: %.6f", model_summary$r_squared)
    )
    writeLines(lines, output_path)
}


write_final_pool_effect_outputs <- function(
    analysis,
    input_data,
    output_dir
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("This endpoint requires ggplot2 for its two figures.", call. = FALSE)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    write.table(
        analysis$pool_effects,
        file.path(output_dir, "pool-effects.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    write.table(
        analysis$model_summary,
        file.path(output_dir, "model-summary.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    write_final_model_summary(
        input_data,
        analysis$model_summary,
        file.path(output_dir, "analysis-summary.txt")
    )
    ggplot2::ggsave(
        file.path(output_dir, "pool-effects.png"),
        make_final_pool_effect_plot(analysis$pool_effects),
        width = 9,
        height = 7,
        dpi = 300
    )
    ggplot2::ggsave(
        file.path(output_dir, "bird-vs-birdbath-by-pool.png"),
        make_final_scatter_plot(analysis, input_data$pool_levels),
        width = 13,
        height = 10,
        dpi = 300
    )
}
