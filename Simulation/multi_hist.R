library(dplyr)
library(ggplot2)
library(tidyr)


multi_hist <- function(..., Time_list, theta_coef,
                       dim_obs_list = NULL,
                       source_names = NULL,
                       binwidth = 0.02,
                       compare_p = FALSE, # compare wrt p (dim_obs)
                       compare_method = NULL, # method to compare different p
                       compare_T = NULL, # fix T to compare different p
                       compare_dim_obs = NULL, # choose which dim combos to compare
                       p_ncol = 2, # layout in facet_wrap
                       rep = 100,
                       acc = 2, # accuracy coef (blue dotted lines)
                       compare_dataset = FALSE, # compare whole vs extracted
                       dataset_names = NULL, # e.g. c("All","Accurate") in the same order as inputs
                       dataset_alpha = 0.55, # transparency when overlaying
                       dataset_position = c("identity","dodge") # how datasets are arranged
) {
  
  dataset_position <- match.arg(dataset_position)
  
  dots <- list(...)
  
  dfs <- lapply(seq_along(dots), function(i) {
    x <- dots[[i]]
    df <- if (!is.null(x$cp_est_df)) x$cp_est_df else x
    if (!is.data.frame(df)) return(NULL)
    
    nm <- if (!is.null(dataset_names) && length(dataset_names) >= i) {
      dataset_names[[i]]
    } else {
      paste0("data", i)
    }
    df$Dataset <- nm
    df
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (!length(dfs)) stop("Pass at least one object containing cp_est_df (or a cp_est_df data.frame).")
  
  result <- dplyr::bind_rows(dfs)
  
  if (!is.null(source_names)) {
    result <- dplyr::filter(result, method %in% source_names)
  }
  if (!is.null(dim_obs_list)) {
    dim_keep <- vapply(dim_obs_list, paste, collapse = "x", FUN.VALUE = character(1))
    result <- dplyr::filter(result, dim_obs %in% dim_keep)
  }
  
  time_labels <- paste0("T=", Time_list)
  
  make_true_cp <- function(Tvals) {
    tibble::tibble(Tval = Tvals) %>%
      dplyr::mutate(
        w  = round(acc * log(Tval)),
        i1 = floor(theta_coef[1] * Tval),
        i2 = floor(theta_coef[2] * Tval),
        i3 = floor(theta_coef[3] * Tval),
        sc1 = i1 / Tval,
        sc2 = i2 / Tval,
        sc3 = i3 / Tval,
        l1 = pmax(0, (i1 - w) / Tval),
        r1 = pmin(1, (i1 + w) / Tval),
        l2 = pmax(0, (i2 - w) / Tval),
        r2 = pmin(1, (i2 + w) / Tval),
        l3 = pmax(0, (i3 - w) / Tval),
        r3 = pmin(1, (i3 + w) / Tval),
        Time = factor(Tval, levels = Time_list, labels = time_labels)
      ) %>%
      dplyr::select(Time, dplyr::starts_with("sc"), dplyr::starts_with("l"), dplyr::starts_with("r")) %>%
      tidyr::pivot_longer(cols = -Time, names_to = "key", values_to = "x") %>%
      dplyr::mutate(
        cp_name = gsub("^[slr]", "sc", key),
        line_type = dplyr::case_when(
          substr(key, 1, 2) == "sc" ~ "centre",
          substr(key, 1, 1) == "l" ~ "left",
          TRUE ~ "right"
        )
      )
  }
  
  do_compare_dataset <- isTRUE(compare_dataset) && dplyr::n_distinct(result$Dataset) > 1L
  
  # ==========================================================
  # Branch 1: common multiple figs
  # ==========================================================
  if (!isTRUE(compare_p)) {
    
    cp_est_df_plot <- result %>%
      dplyr::filter(!is.na(scaled_cp), Time %in% Time_list) %>%
      dplyr::mutate(
        Time = factor(Time, levels = Time_list, labels = time_labels),
        Method = factor(method, levels = unique(method)),
        Dataset = factor(Dataset, levels = unique(Dataset))
      )
    
    true_cp_df <- make_true_cp(Time_list)
    true_cp_centre <- dplyr::filter(true_cp_df, line_type == "centre")
    true_cp_bounds <- dplyr::filter(true_cp_df, line_type != "centre")
    
    p <- ggplot2::ggplot(cp_est_df_plot, ggplot2::aes(x = scaled_cp))
    
    if (do_compare_dataset) {
      p <- p +
        ggplot2::geom_histogram(
          ggplot2::aes(y = after_stat(count), fill = Dataset),
          binwidth = binwidth,
          colour = "black",
          position = dataset_position,
          alpha = dataset_alpha
        )
    } else {
      p <- p +
        ggplot2::geom_histogram(
          ggplot2::aes(y = after_stat(count)),
          binwidth = binwidth,
          colour = "black",
          position = "identity",
          alpha = 0.7
        )
    }
    
    x_lab <- if (!is.null(dim_obs_list) && length(dim_obs_list) == 1) {
      paste0(
        "Estimated change points, (p1,p2,p3) = (",
        paste(dim_obs_list[[1]], collapse = ","),
        ")"
      )
    } else {
      "Estimated change points"
    }
    
    p <- p +
      ggplot2::geom_vline(
        data = true_cp_centre,
        ggplot2::aes(xintercept = x, group = cp_name),
        linetype = "dashed", colour = "red", linewidth = 0.6
      ) +
      ggplot2::geom_vline(
        data = true_cp_bounds,
        ggplot2::aes(xintercept = x, group = interaction(cp_name, line_type)),
        linetype = "dotted", colour = "steelblue1", linewidth = 0.4
      ) +
      ggplot2::facet_grid(rows = ggplot2::vars(Method), cols = ggplot2::vars(Time), drop = TRUE) +
      ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      ggplot2::scale_y_continuous(limits = c(0, rep), breaks = pretty(c(0, rep)), oob = scales::squish) +
      ggplot2::labs(x = x_lab, y = "Frequencies") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black"),
        panel.grid.major = ggplot2::element_line(colour = "grey90"),
        panel.grid.minor = ggplot2::element_line(colour = "grey95"),
        strip.background = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_text(face = "bold", hjust = 0.5),
        strip.text.y = ggplot2::element_text(face = "bold"),
        legend.position = if (do_compare_dataset) "bottom" else "none"
      )
    
    return(p)
  }
  
  # ==========================================================
  # Branch 2: compare p only, fix T
  # ==========================================================
  if (is.null(compare_method)) stop("compare_p=TRUE: please provide compare_method, e.g. 'TFMseg' or 'LR'.")
  if (is.null(compare_T)) stop("compare_p=TRUE: please provide compare_T, e.g. 800.")
  
  if (!is.null(compare_dim_obs)) {
    dim_keep <- vapply(compare_dim_obs, paste, collapse = "x", FUN.VALUE = character(1))
  } else if (!is.null(dim_obs_list)) {
    dim_keep <- vapply(dim_obs_list, paste, collapse = "x", FUN.VALUE = character(1))
  } else {
    dim_keep <- sort(unique(result$dim_obs))
  }
  
  cp_est_df_plot <- result %>%
    dplyr::filter(
      !is.na(scaled_cp),
      method == compare_method,
      Time == compare_T,
      dim_obs %in% dim_keep
    ) %>%
    dplyr::mutate(
      Time = factor(Time, levels = Time_list, labels = time_labels),
      Method = factor(method, levels = compare_method),
      dim_obs = factor(dim_obs, levels = dim_keep),
      Dataset = factor(Dataset, levels = unique(Dataset))
    )
  
  true_cp_df <- make_true_cp(compare_T)
  true_cp_centre <- dplyr::filter(true_cp_df, line_type == "centre")
  true_cp_bounds <- dplyr::filter(true_cp_df, line_type != "centre")
  
  p <- ggplot2::ggplot(cp_est_df_plot, ggplot2::aes(x = scaled_cp))
  
  if (do_compare_dataset) {
    p <- p +
      ggplot2::geom_histogram(
        ggplot2::aes(y = after_stat(count), fill = Dataset),
        binwidth = binwidth,
        colour = "black",
        position = dataset_position,
        alpha = dataset_alpha
      )
  } else {
    p <- p +
      ggplot2::geom_histogram(
        ggplot2::aes(y = after_stat(count)),
        binwidth = binwidth,
        colour = "black",
        position = "identity",
        alpha = 0.7
      )
  }
  
  p +
    ggplot2::geom_vline(
      data = true_cp_centre,
      ggplot2::aes(xintercept = x, group = cp_name),
      linetype = "dashed", colour = "red", linewidth = 0.6
    ) +
    ggplot2::geom_vline(
      data = true_cp_bounds,
      ggplot2::aes(xintercept = x, group = interaction(cp_name, line_type)),
      linetype = "dotted", colour = "steelblue1", linewidth = 0.4
    ) +
    ggplot2::facet_wrap(~ dim_obs, ncol = p_ncol, drop = TRUE) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, rep), breaks = pretty(c(0, rep)), oob = scales::squish) +
    ggplot2::labs(
      x = "Estimated change points (scaled)",
      y = "Repetitions",
      title = paste0("Method: ", compare_method, ", T=", compare_T)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      panel.grid.minor = ggplot2::element_line(colour = "grey95"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = if (do_compare_dataset) "bottom" else "none"
    )
}