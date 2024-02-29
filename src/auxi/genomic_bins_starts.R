genomic_bins_starts <- function(start = 1, end = 0, bin_number = 0, bin_size = 0) {
  if (bin_number > 0 && bin_size > 0) stop("genomic bins starts: Use either bin number or bin size")
  if (bin_number == 0 && bin_size == 0) stop("genomic bins starts: Use either bin number or bin size")
  if (end < start) stop("genomic bins starts: End smaller than start will not work too well...")
  if (bin_size >= (end - start)) return(start)

  if (bin_number > 0) {
    seq_per_bin <- (end - start + 1) %/% bin_number
    remaining_seq <- (end - start + 1) %% bin_number
    bin_sizes <- rep(seq_per_bin, bin_number)
    if (remaining_seq > 0) {
      add_remaining_here <- sample(1:bin_number, remaining_seq)
      bin_sizes[add_remaining_here] <- bin_sizes[add_remaining_here] + 1
    }
    start_positions <- bin_sizes
    for (i in 2 : length(start_positions)) {
      start_positions[i] <- start_positions[i] + start_positions[i - 1]
    }
    start_positions <- start_positions - start_positions[1] + start
    remove(seq_per_bin, remaining_seq)

    return(start_positions)
  }
  if (bin_size > 0) {
    if ((end - start) < bin_size) return(start)
    start_positions <- seq(start, (end - bin_size), bin_size)
    if ((end - start_positions[length(start_positions)]) < (bin_size / 2)) { #if future last win length is less than half a bin size
      start_positions <- start_positions[-length(start_positions)] # remove the last one
    }
    return(start_positions)
  }
  return(NA)
}
