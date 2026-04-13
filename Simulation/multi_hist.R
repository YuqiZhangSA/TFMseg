library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(scales)


multi_hist <- function(..., Time_list, theta_coef,
                       dim_obs_list = NULL,
                       source_names = NULL,
                       binwidth = 0.0625,
                       compare_p = FALSE, # compare wrt p (dim_obs)
                       compare_method = NULL, # method to compare different p
                       compare_T = NULL, # fix T to compare different p
                       compare_dim_obs = NULL, # choose which dim combos to compare
                       p_ncol = 2, # layout in facet_wrap
                       rep = 100,
                       acc = 2, # accuracy coef (blue dotted lines)
                       compare_dataset = FALSE, # compare whole vs extracted
                       dataset_names = c("Complete","Extract"), # e.g. c("Complete","Extract") in the same order as inputs
                       dataset_alpha = 0.5, # transparency when overlaying
                       dataset_position = "identity" #c("identity","dodge") # how datasets are arranged
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
          center = theta_coef[1],
          colour = "black",
          position = dataset_position,
          alpha = dataset_alpha
        )
    } else {
      p <- p +
        ggplot2::geom_histogram(
          ggplot2::aes(y = after_stat(count)),
          binwidth = binwidth,
          center = theta_coef[1],
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



multi_hist_subplot <- function(..., Time_list, theta_coef,
                               dim_obs_list = NULL,
                               source_names = c("TFMseg", "LR"),
                               binwidth = 0.0625,
                               rep = 100,
                               acc = 2,
                               compare_dataset = TRUE,
                               dataset_names = c("Complete", "Extract"),
                               dataset_alpha = 0.5,
                               dataset_position = c("identity", "dodge"),
                               setting = c("S1", "S2", "none")) {
  
  dataset_position <- match.arg(dataset_position)
  setting <- match.arg(setting)
  
  dots <- list(...)
  
  get_df <- function(x) {
    if (is.data.frame(x)) {
      return(x)
    }
    if (is.list(x) && !is.null(x$cp_est_df) && is.data.frame(x$cp_est_df)) {
      return(x$cp_est_df)
    }
    return(NULL)
  }
  
  dfs <- lapply(seq_along(dots), function(i) {
    x <- dots[[i]]
    df <- get_df(x)
    if (is.null(df)) return(NULL)
    
    nm <- if (!is.null(dataset_names) && length(dataset_names) >= i) {
      dataset_names[[i]]
    } else {
      paste0("data", i)
    }
    df$Dataset <- nm
    df
  })
  
  dfs <- Filter(Negate(is.null), dfs)
  if (!length(dfs)) {
    stop("Pass at least one data.frame or one object containing a data.frame named cp_est_df.")
  }
  
  result <- dplyr::bind_rows(dfs)
  
  if (!is.null(source_names)) {
    result <- dplyr::filter(result, method %in% source_names)
  }
  
  if (!is.null(dim_obs_list)) {
    dim_keep <- vapply(dim_obs_list, paste, collapse = "x", FUN.VALUE = character(1))
    result <- dplyr::filter(result, dim_obs %in% dim_keep)
  } else {
    dim_keep <- sort(unique(result$dim_obs))
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
  
  cp_est_df_plot <- result %>%
    dplyr::filter(!is.na(scaled_cp), Time %in% Time_list) %>%
    dplyr::mutate(
      Time = factor(Time, levels = Time_list, labels = time_labels),
      Method = factor(method, levels = source_names),
      dim_obs = factor(dim_obs, levels = dim_keep),
      Dataset = factor(Dataset, levels = unique(Dataset))
    )
  
  true_cp_df <- make_true_cp(Time_list)
  true_cp_centre <- dplyr::filter(true_cp_df, line_type == "centre")
  true_cp_bounds <- dplyr::filter(true_cp_df, line_type != "centre")
  
  dim_labels <- setNames(
    paste0("(", gsub("x", ",", dim_keep), ")"),
    dim_keep
  )
  
  if (setting == "none") {
    x_lab <- paste0("Estimated change points")
  } else {
    x_lab <- paste0("Estimated change points (", setting, ")")
  }
  
  p <- ggplot2::ggplot(cp_est_df_plot, ggplot2::aes(x = scaled_cp))
  
  if (do_compare_dataset) {
    p <- p +
      ggplot2::geom_histogram(
        ggplot2::aes(y = after_stat(count), fill = Dataset),
        binwidth = binwidth,
        center = theta_coef[1],
        colour = "black",
        position = dataset_position,
        alpha = dataset_alpha
      )
  } else {
    p <- p +
      ggplot2::geom_histogram(
        ggplot2::aes(y = after_stat(count)),
        binwidth = binwidth,
        center = theta_coef[1],
        colour = "black",
        position = "identity",
        alpha = 0.7
      )
  }
  
  p +
    ggplot2::geom_vline(
      data = true_cp_centre,
      ggplot2::aes(xintercept = x, group = cp_name),
      linetype = "dashed",
      colour = "red",
      linewidth = 0.6,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_vline(
      data = true_cp_bounds,
      ggplot2::aes(xintercept = x, group = interaction(cp_name, line_type)),
      linetype = "dotted",
      colour = "steelblue1",
      linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(dim_obs),
      cols = ggplot2::vars(Time, Method),
      labeller = ggplot2::labeller(dim_obs = dim_labels)
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, rep),
      breaks = pretty(c(0, rep)),
      oob = scales::squish
    ) +
    ggplot2::labs(
      x = x_lab,
      y = "Frequencies"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      panel.grid.minor = ggplot2::element_line(colour = "grey95"),
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(face = "bold"),
      strip.text.y = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 14)),
      legend.position = if (do_compare_dataset) "bottom" else "none"
    )
}


multi_hist_subplot_method <- function(
    complete_1,
    extract_1,
    complete_2,
    extract_2,
    Time_list,
    theta_coef,
    dim_obs_list = NULL,
    source_names = c("TFMseg"),
    binwidth = 0.0625,
    rep = 100,
    acc = 2,
    scenario_names = c("Complete data", "Missingness"),
    dataset_names = c("Complete", "Extract"),
    dataset_alpha = 0.5,
    dataset_position = "identity",
    setting = c("S1", "S2")
) {
  
  dataset_position <- match.arg(dataset_position, c("identity", "dodge"))
  setting <- match.arg(setting)
  
  get_df <- function(x) {
    if (!is.null(x$cp_est_df)) x$cp_est_df else x
  }
  
  df1 <- get_df(complete_1)
  df2 <- get_df(extract_1)
  df3 <- get_df(complete_2)
  df4 <- get_df(extract_2)
  
  if (!all(vapply(list(df1, df2, df3, df4), is.data.frame, logical(1)))) {
    stop("All inputs must be data.frames or objects containing cp_est_df.")
  }
  
  df1 <- df1 %>% mutate(Scenario = scenario_names[1], Dataset = dataset_names[1])
  df2 <- df2 %>% mutate(Scenario = scenario_names[1], Dataset = dataset_names[2])
  df3 <- df3 %>% mutate(Scenario = scenario_names[2], Dataset = dataset_names[1])
  df4 <- df4 %>% mutate(Scenario = scenario_names[2], Dataset = dataset_names[2])
  
  result <- bind_rows(df1, df2, df3, df4)
  
  if (!is.null(source_names)) {
    result <- result %>% filter(method %in% source_names)
  }
  
  if (!is.null(dim_obs_list)) {
    dim_keep <- vapply(dim_obs_list, paste, collapse = "x", FUN.VALUE = character(1))
    result <- result %>% filter(dim_obs %in% dim_keep)
  } else {
    dim_keep <- sort(unique(result$dim_obs))
  }
  
  time_labels <- paste0("T=", Time_list)
  
  make_true_cp <- function(Tvals) {
    tibble(Tval = Tvals) %>%
      mutate(
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
      select(Time, starts_with("sc"), starts_with("l"), starts_with("r")) %>%
      pivot_longer(cols = -Time, names_to = "key", values_to = "x") %>%
      mutate(
        cp_name = gsub("^[slr]", "sc", key),
        line_type = case_when(
          substr(key, 1, 2) == "sc" ~ "centre",
          substr(key, 1, 1) == "l" ~ "left",
          TRUE ~ "right"
        )
      )
  }
  
  cp_est_df_plot <- result %>%
    filter(!is.na(scaled_cp), Time %in% Time_list) %>%
    mutate(
      Time = factor(Time, levels = Time_list, labels = time_labels),
      Method = factor(method, levels = source_names),
      dim_obs = factor(dim_obs, levels = dim_keep),
      Dataset = factor(Dataset, levels = dataset_names),
      Scenario = factor(Scenario, levels = scenario_names)
    )
  
  true_cp_df <- make_true_cp(Time_list)
  true_cp_centre <- true_cp_df %>% filter(line_type == "centre")
  true_cp_bounds <- true_cp_df %>% filter(line_type != "centre")
  
  dim_labels <- setNames(
    paste0("(", gsub("x", ",", dim_keep), ")"),
    dim_keep
  )
  
  x_lab <- paste0("Estimated change points (", setting, ")")
  
  p <- ggplot(cp_est_df_plot, aes(x = scaled_cp))
  
  p <- p +
    geom_histogram(
      aes(y = after_stat(count), fill = Dataset),
      binwidth = binwidth,
      center = theta_coef[1],
      colour = "black",
      position = dataset_position,
      alpha = dataset_alpha
    ) +
    geom_vline(
      data = true_cp_centre,
      aes(xintercept = x, group = cp_name),
      linetype = "dashed",
      colour = "red",
      linewidth = 0.6,
      inherit.aes = FALSE
    ) +
    geom_vline(
      data = true_cp_bounds,
      aes(xintercept = x, group = interaction(cp_name, line_type)),
      linetype = "dotted",
      colour = "steelblue1",
      linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    facet_grid(
      rows = vars(dim_obs),
      cols = vars(Time, Scenario),
      labeller = labeller(dim_obs = dim_labels)
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ) +
    scale_y_continuous(
      limits = c(0, rep),
      breaks = pretty(c(0, rep)),
      oob = squish
    ) +
    labs(
      x = x_lab,
      y = "Frequencies"
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey95"),
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold"),
      strip.text.y = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 14)),
      legend.position = "bottom"
    )
  
  return(p)
}



plot_runtime <- function(res_obj_list,
                              source_labels = NULL,
                              Time_list = c(400, 800, 1600, 3200),
                              dim_obs_list = list(c(10, 10, 10),
                                                  c(10, 10, 100),
                                                  c(10, 20, 40),
                                                  c(20, 20, 20)),
                              runtime_col = "Mean_time_s",
                              error_col = "SD_time_s",
                              raw_runtime_col = "time",
                              method_order = c("TFMseg", "LR"),
                              method_labels = c("TFMseg", "LR"),
                              y_lab = "Runtime (s)",
                              log_y = TRUE,
                              error_style = c("errorbar", "dotted", "ribbon", "none"),
                              point_size = 1.5,
                              line_width = 0.5,
                              errorbar_width = 0.10,
                              dotted_line_width = 0.6,
                              ribbon_alpha = 0.3,
                              base_size = 11,
                              legend_title = NULL,
                              x_text_angle = 25,
                              x_text_size = 11) {
  
  error_style <- match.arg(error_style)
  
  library(ggplot2)
  library(dplyr)
  
  if (!is.list(res_obj_list) || length(res_obj_list) == 0L) {
    stop("`res_obj_list` must be a non-empty list.")
  }
  
  if (is.null(source_labels)) {
    source_labels <- paste0("Setting ", seq_along(res_obj_list))
  }
  
  if (length(source_labels) != length(res_obj_list)) {
    stop("`source_labels` must have the same length as `res_obj_list`.")
  }
  
  make_dim_chr <- function(dim_obs_list) {
    vapply(dim_obs_list, function(x) paste(x, collapse = "x"), character(1))
  }
  
  make_dim_lab <- function(dim_obs_list) {
    vapply(dim_obs_list, function(x) paste0("(", paste(x, collapse = ","), ")"), character(1))
  }
  
  build_runtime_summary <- function(res_obj,
                                    source_label,
                                    runtime_col,
                                    error_col,
                                    raw_runtime_col) {
    
    if (!is.null(res_obj$summary)) {
      df_sum <- as.data.frame(res_obj$summary)
    } else {
      df_sum <- NULL
    }
    
    needed_sum <- c("Time", "dim_obs", "method", runtime_col)
    has_summary_cols <- !is.null(df_sum) && all(needed_sum %in% names(df_sum))
    
    if (has_summary_cols) {
      out <- df_sum
    } else {
      if (is.null(res_obj$df_raw)) {
        stop("Neither usable `summary` nor `df_raw` found in one input object.")
      }
      df_raw <- as.data.frame(res_obj$df_raw)
      req_raw <- c("Time", "dim_obs", "method", raw_runtime_col)
      miss_raw <- setdiff(req_raw, names(df_raw))
      if (length(miss_raw) > 0L) {
        stop("`df_raw` is missing: ", paste(miss_raw, collapse = ", "))
      }
      
      out <- df_raw %>%
        group_by(Time, dim_obs, method) %>%
        summarise(
          !!runtime_col := mean(.data[[raw_runtime_col]], na.rm = TRUE),
          !!error_col   := stats::sd(.data[[raw_runtime_col]], na.rm = TRUE),
          .groups = "drop"
        )
    }
    
    if (!error_col %in% names(out)) {
      out[[error_col]] <- NA_real_
    }
    
    out$source_row <- source_label
    out
  }
  
  dim_obs_chr <- make_dim_chr(dim_obs_list)
  dim_obs_lab <- make_dim_lab(dim_obs_list)
  
  df_list <- Map(
    f = function(obj, lab) {
      build_runtime_summary(
        res_obj = obj,
        source_label = lab,
        runtime_col = runtime_col,
        error_col = error_col,
        raw_runtime_col = raw_runtime_col
      )
    },
    obj = res_obj_list,
    lab = source_labels
  )
  
  df_all <- bind_rows(df_list)
  
  df_all <- df_all %>%
    filter(Time %in% Time_list, dim_obs %in% dim_obs_chr)
  
  if (nrow(df_all) == 0L) {
    stop("No rows remain after filtering by `Time_list` and `dim_obs_list`.")
  }
  
  df_all <- df_all %>%
    mutate(
      Time = factor(Time, levels = Time_list, labels = as.character(Time_list)),
      dim_obs = factor(dim_obs, levels = dim_obs_chr, labels = dim_obs_lab),
      method = factor(method, levels = method_order, labels = method_labels),
      source_row = factor(source_row, levels = source_labels)
    )
  
  has_error <- all(c(runtime_col, error_col) %in% names(df_all)) &&
    any(!is.na(df_all[[error_col]]))
  
  if (has_error) {
    df_all <- df_all %>%
      mutate(
        ymin_runtime = pmax(.data[[runtime_col]] - .data[[error_col]], 1e-8),
        ymax_runtime = .data[[runtime_col]] + .data[[error_col]]
      )
  } else {
    error_style <- "none"
  }
  
  pd <- position_dodge(width = 0.18)
  
  n_source <- length(source_labels)
  use_row_facet <- n_source > 1L
  
  p <- ggplot(
    df_all,
    aes(
      x = Time,
      y = .data[[runtime_col]],
      colour = method,
      shape = method,
      linetype = method,
      group = method
    )
  )
  
  if (error_style == "ribbon" && has_error) {
    ribbon_df <- df_all %>%
      mutate(Time_id = match(Time, levels(Time)))
    
    p <- p +
      geom_ribbon(
        data = ribbon_df,
        aes(
          x = Time_id,
          ymin = ymin_runtime,
          ymax = ymax_runtime,
          fill = method,
          group = method
        ),
        inherit.aes = FALSE,
        alpha = ribbon_alpha,
        colour = NA
      )
  }
  
  p <- p +
    geom_line(linewidth = line_width, position = pd) +
    geom_point(size = point_size, position = pd)
  
  if (error_style == "errorbar" && has_error) {
    p <- p +
      geom_errorbar(
        aes(ymin = ymin_runtime, ymax = ymax_runtime),
        width = errorbar_width,
        position = pd
      )
  }
  
  if (error_style == "dotted" && has_error) {
    p <- p +
      geom_line(
        aes(y = ymin_runtime),
        linewidth = dotted_line_width,
        alpha = 0.6,
        position = pd,
        show.legend = FALSE
      ) +
      geom_line(
        aes(y = ymax_runtime),
        linewidth = dotted_line_width,
        alpha = 0.6,
        position = pd,
        show.legend = FALSE
      )
  }
  
  p <- p +
    geom_hline(
      yintercept = c(1, 10, 100),
      linewidth = 0.25,
      linetype = "dashed",
      alpha = 0.6
    )
  
  if (use_row_facet) {
    p <- p +
      facet_grid(
        rows = vars(source_row),
        cols = vars(dim_obs),
        labeller = labeller(
          source_row = label_parsed,
          dim_obs = label_value
        )
      )
  } else {
    p <- p +
      facet_grid(
        cols = vars(dim_obs),
        labeller = labeller(dim_obs = label_value)
      )
  }
  
  p <- p +
    labs(
      x = "Time",
      y = y_lab,
      colour = legend_title,
      shape = legend_title,
      linetype = legend_title,
      fill = legend_title
    ) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = x_text_angle,
        hjust = 1,
        vjust = 1,
        size = x_text_size
      )
    ) +
    scale_shape_manual(values = c(16, 17)) +
    scale_linetype_manual(values = c("solid", "dashed"))
  
  if (error_style == "ribbon" && has_error) {
    p <- p +
      scale_x_discrete(drop = FALSE) +
      scale_fill_discrete(guide = "none")
  }
  
  if (log_y) {
    p <- p + scale_y_log10()
  }
  
  p
}