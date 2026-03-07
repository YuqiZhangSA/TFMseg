seeded_intervals <- function(Time, minl = NULL) {
  if (Time < 3) stop("Time should be at least 3")
  Time <- as.integer(Time)
  Time_d <- as.double(Time)
  
  if (is.null(minl)) {
    minl <- as.integer(0.5 * floor(Time_d / log(Time_d)))
  } else {
    minl <- as.integer(minl)
  }
  
  # choose layers 2*m_h >= Time/log(Time)
  h_max <- as.integer(floor(log(4 * log(Time_d), base = 2)))
  if (is.na(h_max) || h_max < 1L) {
    return(data.frame(st = integer(0), ed = integer(0)))
  }
  
  out_st <- integer(0)
  out_ed <- integer(0)
  
  for (h in 1:h_max) {
    m_h <- Time_d / (2^h)  
    
    i_max <- as.integer(floor(Time_d / m_h) - 1L)
    if (i_max < 1L) next
    i <- seq_len(i_max)
    
    st <- floor((i - 1) * m_h)
    ed <- ceiling((i + 1) * m_h)
    
    st[st < 0] <- 0
    ed[ed > Time] <- Time
    
    keep <- (ed - st) >= minl
    out_st <- c(out_st, st[keep])
    out_ed <- c(out_ed, ed[keep])
  }
  
  intervals <- unique(data.frame(st = as.integer(out_st), ed = as.integer(out_ed)))
  intervals <- intervals[order(intervals$st, intervals$ed), , drop = FALSE]
  rownames(intervals) <- NULL
  intervals
}



# for (Time in c(400, 800, 1600, 3200)) {
#   iv <- seeded_intervals(Time, floor(0.5*Time/log(Time)))
#   cat("Time=", Time,
#       " #intervals=", nrow(iv),
#       " minLen=", min(iv$ed - iv$st),
#       " 0.5*Time/log(Time)=", floor(0.5 * Time/log(Time)), "\n")
# }


