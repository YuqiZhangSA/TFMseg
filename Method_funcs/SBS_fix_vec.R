seeded_intervals_vec <- function(Time, minl = 2L) {
  if (Time < 3) stop("T should be at least 3")
  Time <- as.integer(Time)
  minl <- as.integer(minl)
  
  mu_T <- as.integer(floor(log2(Time)))  
  
  out <- vector("list", mu_T)
  
  for (h in seq_len(mu_T)) {
    m_h <- Time / (2^h)                  
    i_max <- as.integer(floor(Time / m_h) - 1L)
    
    if (i_max < 1L) {
      out[[h]] <- NULL
      next
    }
    
    i <- seq_len(i_max)
    
    st <- floor((i - 1) * m_h) + 1L   
    ed <- ceiling((i + 1) * m_h)      
    
    st[st < 1L] <- 1L
    ed[ed > Time]  <- Time
    
    keep <- (ed - st + 1L) >= minl
    out[[h]] <- data.frame(st = st[keep], ed = ed[keep])
  }
  
  intervals <- do.call(rbind, out)
  intervals <- unique(intervals)
  intervals <- intervals[order(intervals$st, intervals$ed), , drop = FALSE]
  rownames(intervals) <- NULL
  intervals
}

# for (Time in c(400, 800, 1600, 3200)) {
#   iv <- seeded_intervals(Time, floor(Time/log(Time)))
#   cat("T=", Time,
#       " #intervals=", nrow(iv),
#       " minLen=", min(iv$ed - iv$st),
#       " T/log(T)=", floor(Time/log(Time)), "\n")
# }
