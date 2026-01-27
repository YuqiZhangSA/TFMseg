library(dplyr)
library(ggplot2)
library(tidyr)

single_hist_clean <- function(..., Time_list, theta_coef) {
  dots <- list(...)
  dfs  <- lapply(dots, function(x) if (!is.null(x$cp_est_df)) x$cp_est_df else x)
  stopifnot(all(vapply(dfs, is.data.frame, TRUE)))
  result <- bind_rows(dfs)
  
  time_labels <- paste0("T=", Time_list)
  
  cp_est_df_plot <- result %>%
    filter(!is.na(scaled_cp), Time %in% Time_list) %>%
    mutate(
      Time   = factor(Time, levels = Time_list, labels = time_labels),
      Method = case_when(
        Group == "LR" ~ "LR",
        Group == "TFMseg: fixed, diag" ~ "TFMseg: fixed, diag",
        Group == "TFMseg: fixed, full" ~ "TFMseg: fixed, full",
        Group == "FMseg: oracle, diag" ~ "FMseg: oracle, diag",
        Group == "FMseg: fixed, diag"  ~ "FMseg: fixed, diag",
        Group == "FMseg: oracle, full" ~ "FMseg: oracle, full",
        Group == "FMseg: fixed, full"  ~ "FMseg: fixed, full",
        TRUE ~ as.character(Group)
      )
    )
  
  preferred_levels <- c("TFMseg: fixed, diag", "TFMseg: fixed, full",
                        "FMseg: oracle, diag", "FMseg: fixed, diag",
                        "FMseg: oracle, full", "FMseg: fixed, full",
                        "LR")
  present_levels <- intersect(preferred_levels, unique(cp_est_df_plot$Method))
  cp_est_df_plot <- cp_est_df_plot %>%
    mutate(Method = factor(Method, levels = present_levels)) %>%
    droplevels()
  
  true_cp_df <- tibble(Tval = Time_list) %>%
    mutate(
      cp1 = floor(theta_coef[1] * Tval),
      cp2 = floor(theta_coef[2] * Tval),
      cp3 = floor(theta_coef[3] * Tval),
      sc1 = cp1 / Tval, sc2 = cp2 / Tval, sc3 = cp3 / Tval,
      Time = factor(Tval, levels = Time_list, labels = time_labels)
    ) %>%
    pivot_longer(starts_with("sc"), names_to = "cp_name", values_to = "scaled_cp")
  
  ggplot(cp_est_df_plot, aes(x = scaled_cp, fill = Method)) +
    geom_histogram(aes(y = ..count.. / 3), binwidth = 0.02,
                   colour = "black", position = "identity", alpha = 0.7) +
    geom_vline(data = true_cp_df, aes(xintercept = scaled_cp),
               linetype = "dotted", colour = "red", size = 0.8) +
    facet_grid(rows = vars(Method), cols = vars(Time), drop = TRUE) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_fill_brewer(palette = "Set2", drop = TRUE) +
    labs(x = "Estimated Change Points", y = "Count") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey95"),
      strip.background = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text.x = element_text(face = "bold", hjust = 0.5),
      strip.text.y = element_text(face = "bold"),
      legend.position = "bottom"
    )
}




multi_hist <- function(..., Time_list, theta_coef,
                       dim_obs_list = NULL,
                       source_names = NULL,
                       binwidth = 0.02) {
  dots <- list(...)
  dfs  <- lapply(dots, function(x) if (!is.null(x$cp_est_df)) x$cp_est_df else x)
  dfs  <- Filter(is.data.frame, dfs)
  if (!length(dfs)) stop("Pass at least one object containing cp_est_df (or a cp_est_df data.frame).")
  
  result <- dplyr::bind_rows(dfs)
  
  # optional filters
  if (!is.null(source_names)) {
    result <- dplyr::filter(result, method %in% source_names)
  }
  if (!is.null(dim_obs_list)) {
    dim_keep <- vapply(dim_obs_list, paste, collapse = "x", FUN.VALUE = character(1))
    result <- dplyr::filter(result, dim_obs %in% dim_keep)
  }
  
  time_labels <- paste0("T=", Time_list)
  
  cp_est_df_plot <- result %>%
    dplyr::filter(!is.na(scaled_cp), Time %in% Time_list) %>%
    dplyr::mutate(
      Time   = factor(Time, levels = Time_list, labels = time_labels),
      Method = factor(method, levels = unique(method))
    )
  
  true_cp_df <- tibble::tibble(Tval = Time_list) %>%
    dplyr::mutate(
      sc1 = floor(theta_coef[1] * Tval) / Tval,
      sc2 = floor(theta_coef[2] * Tval) / Tval,
      sc3 = floor(theta_coef[3] * Tval) / Tval,
      Time = factor(Tval, levels = Time_list, labels = time_labels)
    ) %>%
    tidyr::pivot_longer(dplyr::starts_with("sc"),
                        names_to = "cp_name", values_to = "scaled_cp")
  
  ggplot2::ggplot(cp_est_df_plot, ggplot2::aes(x = scaled_cp)) +
    ggplot2::geom_histogram(ggplot2::aes(y = after_stat(count) / 3),
                            binwidth = binwidth, colour = "black",
                            position = "identity", alpha = 0.7) +
    ggplot2::geom_vline(data = true_cp_df, ggplot2::aes(xintercept = scaled_cp),
                        linetype = "dotted", colour = "red", linewidth = 0.8) +
    ggplot2::facet_grid(rows = ggplot2::vars(Method), cols = ggplot2::vars(Time), drop = TRUE) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::labs(x = "Estimated change points (scaled)", y = "Count") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      panel.grid.minor = ggplot2::element_line(colour = "grey95"),
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.text.y = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
}

