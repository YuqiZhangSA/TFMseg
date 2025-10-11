align_change_points <- function(estimated_points, true_points, Time, tolerance = 50) {
  n_true <- length(true_points)
  result <- rep(NA, n_true)
  used <- rep(FALSE, length(estimated_points))
  
  for (i in seq_len(n_true)) {
    candidates <- which(!used & abs(estimated_points - true_points[i]) <= tolerance)
    if (length(candidates) > 0) {
      deviations <- abs(estimated_points[candidates] - true_points[i])
      best <- candidates[which.min(deviations)]
      result[i] <- estimated_points[best] - true_points[i]
      used[best] <- TRUE
    }
  }
  return(result)
}

compute_stats <- function(Time, true_cp, detected_cp_list, method_name = "") {
  nrep <- length(detected_cp_list)
  m_true <- length(true_cp)
  deviations <- vector("list", m_true)
  accuracy_counts_per_cp <- integer(m_true)
  m_est_diff <- integer(nrep)
  for (sim in seq_len(nrep)) {
    est_cp <- detected_cp_list[[sim]]
    m_est_diff[sim] <- length(est_cp) - m_true
    dev <- align_change_points(est_cp, true_cp, Time, tolerance = 2*log(Time))
    for (j in seq_len(m_true)) {
      deviations[[j]] <- c(deviations[[j]], dev[j])
      if (!is.na(dev[j]) && abs(dev[j]) <= 2 * log(Time)) {
        accuracy_counts_per_cp[j] <- accuracy_counts_per_cp[j] + 1
      }
    }
  }
  freq_m_diff_0  <- sum(m_est_diff == 0) / nrep
  freq_m_diff_p1 <- sum(m_est_diff == 1) / nrep
  freq_m_diff_m1 <- sum(m_est_diff == -1) / nrep
  freq_m_diff_ge2 <- sum(m_est_diff >= 2) / nrep
  freq_m_diff_le2 <- sum(m_est_diff <= -2) / nrep
  
  cat("Method:", method_name, "\n")
  cat("Frequency of m_est - m:\n")
  cat("  Value = 0:", freq_m_diff_0, "\n")
  cat("  Value = +1:", freq_m_diff_p1, "\n")
  cat("  Value = -1:", freq_m_diff_m1, "\n")
  cat("  Value >= 2:", freq_m_diff_ge2, "\n")
  cat("  Value <= -2:", freq_m_diff_le2, "\n\n")
  
  detection_rate <- sapply(deviations, function(x) mean(!is.na(x)))
  accuracy_per_cp <- accuracy_counts_per_cp / nrep
  for (i in seq_len(m_true)) {
    cat("Change Point", i, " (True CP =", true_cp[i], "):\n")
    cat("  Detection Rate:", round(detection_rate[i] * 100, 2), "%\n")
    cat("  Accuracy for CP", i, ":", round(accuracy_per_cp[i] * 100, 2), "%\n\n")
  }
  invisible(deviations)
}



visualise_stats <- function(detected_cp_list, true_cp, Time, method_name = "") {
  all_cps <- unlist(detected_cp_list)
  all_cps <- as.numeric(all_cps)
  all_cps <- all_cps[!is.na(all_cps) & is.finite(all_cps)]
  scaled_cps <- all_cps / Time
  scaled_true_cp <- true_cp / Time
  hist(
    scaled_cps,
    breaks = seq(0, 1, by = 2 * (log(Time) / Time)),
    main   = paste("Distribution of Scaled Estimated Change Points\n", method_name),
    xlab   = "Scaled Estimated Change Points",
    col    = "skyblue", border = "black",
    xlim   = c(0, 1)
  )
  abline(v = scaled_true_cp, col = "red", lty = "dotted", lwd = 2)
  for (cp in scaled_true_cp) {
    abline(v = cp - 2 * log(Time) / Time, col = "orange", lty = "dashed")
    abline(v = cp + 2 * log(Time) / Time, col = "orange", lty = "dashed")
  }
}
