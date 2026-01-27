plot_FMseg_iterations <- function(results, Time, m = 3, trim, threshold = NULL, method = c("oracle", "fixed"), theta_coef = c(0.25,0.5,0.75) ){
  #true_cp <- floor(Time * 1:m / (m + 1))
  true_cp <- c(floor(Time * theta_coef[1]), floor(Time * theta_coef[2]), floor(Time * theta_coef[3]))
  method <- match.arg(method, c("fixed", "oracle", "auto"))
  results <- results[order(results$ed - results$st, decreasing = FALSE), ]
  
  all_results <- results  
  selected_results <- list(all_results)
  selected_thresholds <- NULL
  final_threshold <- NULL
  next_highest_thresholds <- list()  
  
  compute_next_threshold <- function(current_threshold, S, thds, all_results, Time) {
    cand_candidates <- sort(thds[thds < current_threshold], decreasing = TRUE)
    
    for (cand in cand_candidates) {
      candidate_intervals <- all_results[all_results$val == cand, ]
      
      if (nrow(candidate_intervals) > 0) {
        for (j in seq_len(nrow(candidate_intervals))) {
          interval <- candidate_intervals[j, ]
          st <- interval$st
          ed <- interval$ed
          
          # Check if any detected CP is within trim distance of the interval
          # if (all(!between(S, interval$st, interval$ed)))
          if (all(S - st < log(Time) | ed - S < log(Time))) {
            return(cand)
          }
        }
      }
    }
    
    return(NULL)  # No valid next-highest threshold found
  }
  
  
  
  if (method == "fixed") {
    if (is.null(threshold)) {
      stop("Must supply 'threshold' for fixed method.")
    }
    
    selected_results[[1]] <- results
    results <- results[results$val >= threshold + 5e-3, ]
    
    #results <- results[order(results$val, decreasing = TRUE), ]
    
    if (nrow(results) == 0) {
      warning("No intervals exceed the given threshold.")
      return(results)
    }
    
    cid <- 0
    selected_indices <- integer()
    while (cid < nrow(results)) {
      cid <- cid + 1
      selected_indices <- c(selected_indices, cid)
      
      if (cid < nrow(results)) {
        tmp <- results[(cid + 1):nrow(results), ]
        tmp <- tmp[!(tmp$st < results$est.cp[cid] & tmp$ed >= results$est.cp[cid]), ]
        results <- rbind(results[1:cid, ], tmp)
      }
    }
    selected_results[[2]] <- results[selected_indices, ]
    
  } else if (method == "oracle") {
    thds <- sort(results$val, decreasing = TRUE)
    threshold_candidates <- list()
    
    for (i in seq_along(thds)) {
      thd <- thds[i]
      aux <- results[results$val >= thd, ]
      cid <- 0
      selected_indices <- integer()
      while (cid < m && cid < nrow(aux)) {
        cid <- cid + 1
        selected_indices <- c(selected_indices, cid)
        if (cid < nrow(aux)) {
          tmp <- aux[(cid + 1):nrow(aux), ]
          tmp <- tmp[!(tmp$st < aux$est.cp[cid] & tmp$ed >= aux$est.cp[cid]), ]
          aux <- rbind(aux[1:cid, ], tmp)
        }
      }
      threshold_candidates[[i]] <- aux
      if (cid == m) {
        selected_thresholds <- thds[1:i]
        final_threshold <- thds[i]
        S_final <- aux$est.cp
        next_val_final <- compute_next_threshold(final_threshold, S_final, thds, all_results, Time)
        next_highest_thresholds[[i]] <- next_val_final
        break
      } else {
        S_candidate <- aux$est.cp
        next_highest_thresholds[[i]] <- compute_next_threshold(thd, S_candidate, thds, all_results, Time)
      }
    }
    
    selected_results <- c(list(all_results), threshold_candidates[1:length(selected_thresholds)])
  }
  
  plots <- list()
  x_limits <- range(all_results$st - 25, all_results$ed + 35, na.rm = TRUE)
  y_limits <- c(0, max(all_results$val, na.rm = TRUE) + 2.5)
  
  for (i in seq_along(selected_results)) {
    df <- selected_results[[i]] %>% mutate(interval_id = row_number())
    
    df_true <- df %>% 
      rowwise() %>% 
      mutate(true_cp_in_interval = list(true_cp[true_cp >= st & true_cp <= ed])) %>% 
      ungroup() %>% 
      unnest(cols = c(true_cp_in_interval)) %>% 
      rename(true_cp = true_cp_in_interval)
    
    half_height <- ifelse(length(y_limits) > 1, diff(y_limits) * 0.02, 0.1)
    
    p <- ggplot() +
      geom_segment(data = df, aes(x = st, xend = ed, y = val, yend = val),
                   linetype = "dotted", colour = "black") +
      geom_point(data = df, aes(x = est.cp, y = val),
                 colour = "blue", size = 2) +
      geom_segment(data = df_true,
                   aes(x = true_cp, xend = true_cp,
                       y = val - half_height, yend = val + half_height),
                   colour = "red", size = 1) +
      geom_text(data = df, aes(x = est.cp, y = val, label = round(val, 2)),
                vjust = -1, colour = "black", size = 3) +
      labs(x = "Time", y = "CUSUM Value",
           title = paste("Iteration", i - 1, "- Threshold Selection")) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE)
    
    if (method == "fixed") {
      p <- p + geom_hline(yintercept = threshold, linetype = "dashed", colour = "darkorange") +
        annotate("text", x = max(all_results$ed) + 5, y = threshold + 0.5, label = round(threshold, 2), colour = "darkorange")
    }
    
    if (method == "oracle") {
      if (i == 1 && !is.null(final_threshold)) {
        next_val_f <- compute_next_threshold(final_threshold, selected_results[[i]]$est.cp, thds, all_results, Time)
        p <- p + geom_hline(yintercept = final_threshold, linetype = "dashed", colour = "darkorange") +
          annotate("text", x = max(all_results$ed) + 5, y = final_threshold + 0.5, label = round(final_threshold, 2), colour = "darkorange")
        if (!is.null(next_val_f)) {
          p <- p + geom_hline(yintercept = next_val_f, linetype = "dashed", colour = "green") +
            annotate("text", x = max(all_results$ed) + 5, y = next_val_f - 0.5, label = round(next_val_f, 2), colour = "green")
        }
      }
      if (i > 1) {
        current_threshold <- selected_thresholds[i - 1]
        S_current <- selected_results[[i]]$est.cp
        next_val <- compute_next_threshold(current_threshold, S_current, thds, all_results, Time)
        p <- p + 
          geom_hline(yintercept = current_threshold, linetype = "dashed", colour = "darkorange") +
          annotate("text", x = max(all_results$ed) + 5, y = current_threshold + 0.5, label = round(current_threshold, 2), colour = "darkorange")
        if (!is.null(next_val)) {
          p <- p + 
            geom_hline(yintercept = next_val, linetype = "dashed", colour = "green") +
            annotate("text", x = max(all_results$ed) + 5, y = next_val - 0.5, label = round(next_val, 2), colour = "green")
        }
      }
    }
    
    print(p)
  }
}