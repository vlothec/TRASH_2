calc_plot_regs <- function(starts.reps, windows.to.check.if.any.repetas.are.in.and.plot) {
  bins.data <- hist(starts.reps, breaks = ceiling((max(starts.reps) - min(starts.reps)) / windows.to.check.if.any.repetas.are.in.and.plot), plot = FALSE)

  bins.starts <- bins.data$breaks[1:(length(bins.data$breaks) - 1)]
  bins.ends <- bins.data$breaks[2:length(bins.data$breaks)]
  bins.count <- bins.data$counts
  bins.df <- data.frame(bins.starts, bins.ends, bins.count)
  bins.df <- bins.df[bins.df$bins.count != 0, ]
  i <- 1
  while (i < nrow(bins.df)) {
    if (bins.df$bins.ends[i] == bins.df$bins.starts[i + 1]) {
      bins.df$bins.ends[i] <- bins.df$bins.ends[i + 1]
      bins.df$bins.count[i] <- bins.df$bins.count[i] + bins.df$bins.count[i + 1]
      bins.df <- bins.df[-(i + 1),]
      i <- i - 1
    }
    i <- i + 1
  }
  bins.df$bin.size <- bins.df$bins.ends - bins.df$bins.starts

  bins.df$dist.to.next <- 0
  bins.df$correction <- bins.df$bins.starts[1]
  i <- 1
  while (i < nrow(bins.df)) {
    bins.df$dist.to.next[i] <- bins.df$bins.starts[i + 1] - bins.df$bins.ends[i]
    bins.df$correction[i + 1] <- bins.df$correction[i + 1] + sum(bins.df$dist.to.next[1:i])
    i <- i + 1
  }
  bins.df$start.adjusted <- bins.df$bins.starts - bins.df$correction
  bins.df$end.adjusted <- bins.df$bins.ends - bins.df$correction

  return(bins.df)
}