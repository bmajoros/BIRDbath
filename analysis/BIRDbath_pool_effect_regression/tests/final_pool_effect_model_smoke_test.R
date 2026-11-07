command_line <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command_line, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[1L]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "models.R"))
source(file.path(repo_root, "R", "input.R"))
source(file.path(repo_root, "R", "data.R"))
source(file.path(repo_root, "R", "final_pool_effect_model.R"))

set.seed(41)
n_variants <- 800L
pool_levels <- paste0("pool", seq_len(4L))
variant <- rep(paste0("v", seq_len(n_variants)), each = length(pool_levels))
pool <- factor(rep(pool_levels, n_variants), levels = pool_levels)
x <- rep(rnorm(n_variants), each = length(pool_levels))
true_alpha <- -0.12
true_gamma <- 0.35
true_beta <- c(-0.18, -0.04, 0.07, 0.15)
y <- true_alpha + true_gamma * x + true_beta[as.integer(pool)] +
    rnorm(length(x), sd = 0.20)
model_data <- data.frame(
    variant = variant,
    pool = pool,
    log2_theta_birdbath = x,
    log2_theta_bird = y
)

analysis <- fit_final_pool_effect_model(model_data, pool_levels)

stopifnot(
    abs(sum(analysis$pool_effects$beta_log2)) < 1e-10,
    abs(analysis$model_summary$alpha - true_alpha) < 0.03,
    abs(analysis$model_summary$gamma - true_gamma) < 0.03,
    max(abs(analysis$pool_effects$beta_log2 - true_beta)) < 0.03,
    nrow(analysis$pool_effects) == length(pool_levels)
)

ready_wide <- data.frame(
    variant = paste0("v", seq_len(n_variants)),
    theta_birdbath = 2^x[seq(1L, length(x), by = length(pool_levels))],
    pool1 = 2^y[pool == "pool1"],
    pool2 = 2^y[pool == "pool2"],
    pool3 = 2^y[pool == "pool3"],
    pool4 = 2^y[pool == "pool4"],
    check.names = FALSE
)
object_input <- read_pool_effect_data(ready_wide)
stopifnot(
    identical(object_input$input_file, "<in-memory data>"),
    nrow(object_input$model_data) == n_variants * length(pool_levels),
    identical(object_input$pool_levels, pool_levels)
)

ready_file <- tempfile(fileext = ".tsv")
write.table(
    ready_wide,
    ready_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)
file_input <- read_pool_effect_data(ready_file)
unlink(ready_file)
stopifnot(
    nrow(file_input$model_data) == nrow(object_input$model_data),
    identical(file_input$pool_levels, object_input$pool_levels)
)

cat("FINAL_POOL_EFFECT_MODEL_SMOKE_TEST_OK\n")
