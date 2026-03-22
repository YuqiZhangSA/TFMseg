#' @param results_dir folder containing many .rds files (each file is one replicate output)
#' @param dep logical; stored in output as a flag
#' @param out_dir optional folder to save summary table (rds + csv). If NULL, no saving.
#' @param out_name file name to save, e.g. "res_indep_rnull"
#' @param strict if TRUE, stop when required fields missing; else fill with NA when possible
#'
#' @return a list with:
#'   - summary: aggregated table (Time x dim_obs x Method)
#'   - df_raw: row-per-replicate metrics
#'   - cp_est_df: long cp estimates for histograms
#'   - check_counts: replicate counts per (Time, dim_obs, Method)

library(dplyr)
library(tidyr)

`%||%` <- function(x, y) if (is.null(x)) y else x

combine_tasks <- function(results_dir,
                          dep = FALSE,
                          out_dir = NULL,
                          out_name = "res_table",
                          strict = TRUE,
                          build_cp = TRUE,
                          progress_every = 200L,
                          method_filter = NULL) {
  
  results_dir <- path.expand(results_dir)
  if (!dir.exists(results_dir)) stop("results_dir does not exist: ", results_dir)
  
  files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No .rds files found in: ", results_dir)
  
  if (!is.null(method_filter)) {
    files <- files[grepl(paste0("_method-", method_filter, "_"), basename(files), fixed = FALSE)]
  }
  
  if (!length(files)) stop("No matching .rds files found after method filtering.")
  
  df_list <- vector("list", length(files))
  cp_list <- if (build_cp) vector("list", length(files)) else NULL
  
  for (i in seq_along(files)) {
    if (!is.null(progress_every) && progress_every > 0L && i %% progress_every == 0L) {
      message("Reading ", i, "/", length(files), "  (", basename(files[i]), ")")
    }
    
    z <- readRDS(files[i])
    
    if (strict) {
      stopifnot(!is.null(z$res), !is.null(z$Time), !is.null(z$dim_obs),
                !is.null(z$method), !is.null(z$rep_id))
    }
    
    res  <- z$res
    freq <- res$freq_summary %||% rep(NA_real_, 5)
    acc  <- res$accuracy_summary %||% rep(NA_real_, 3)
    tsec <- res$time_sec %||% NA_real_
    
    df_list[[i]] <- tibble(
      Time = z$Time,
      dim_obs = paste(z$dim_obs, collapse = "x"),
      method = as.character(z$method),
      rep_id = z$rep_id,
      le_m2 = freq[1],
      m1 = freq[2],
      m0 = freq[3],
      p1 = freq[4],
      ge_p2 = freq[5],
      acc1 = acc[1] %||% NA_real_,
      acc2 = acc[2] %||% NA_real_,
      acc3 = acc[3] %||% NA_real_,
      time = if (is.numeric(tsec)) tsec[1] else NA_real_,
      dep = dep,
      r_hat = z$r_hat %||% NA_real_
    )
    
    if (build_cp) {
      cps <- z$res$cp_est_scaled[[1]] %||% numeric(0)
      cp_list[[i]] <- if (!length(cps)) {
        tibble(
          scaled_cp = NA_real_,
          sim = z$rep_id,
          Time = z$Time,
          dim_obs = paste(z$dim_obs, collapse = "x"),
          method = as.character(z$method)
        )
      } else {
        tibble(
          scaled_cp = cps,
          sim = z$rep_id,
          Time = z$Time,
          dim_obs = paste(z$dim_obs, collapse = "x"),
          method = as.character(z$method)
        )
      }
    }
  }
  
  df_raw <- bind_rows(df_list)
  cp_est_df <- if (build_cp) bind_rows(cp_list) else NULL
  
  check_counts <- df_raw %>%
    count(Time, dim_obs, method, name = "n") %>%
    arrange(Time, dim_obs, method)
  
  summary <- df_raw %>%
    group_by(Time, dim_obs, method) %>%
    summarise(
      `<= -2` = mean(le_m2, na.rm = TRUE),
      `-1` = mean(m1, na.rm = TRUE),
      `0` = mean(m0, na.rm = TRUE),
      `1` = mean(p1, na.rm = TRUE),
      `>= 2` = mean(ge_p2, na.rm = TRUE),
      Accuracy_j1 = mean(acc1, na.rm = TRUE),
      Accuracy_j2 = mean(acc2, na.rm = TRUE),
      Accuracy_j3 = mean(acc3, na.rm = TRUE),
      Mean_time_s = mean(time, na.rm = TRUE),
      SD_time_s = sd(time, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Time, dim_obs, method)
  
  if (!is.null(out_dir)) {
    out_dir <- path.expand(out_dir)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(summary, file.path(out_dir, paste0(out_name, ".rds")))
    write.csv(summary, file.path(out_dir, paste0(out_name, ".csv")), row.names = FALSE)
  }
  
  list(
    summary = summary,
    df_raw = df_raw,
    cp_est_df = cp_est_df,
    check_counts = check_counts,
    meta = list(results_dir = results_dir, dep = dep, n_files = length(files), build_cp = build_cp)
  )
}
