sequence_window_score <- function(fasta_sequence, window_size, log_messages) {
  # TODO check if changing these settings below can make the script work better,
  # although these were optimised
  kmer <- 12
  fraction_p <- 0.5
  if((window_size / 2) <= kmer) stop("sequence_window_score: window size is too small")

  sequence_full_length <- length(fasta_sequence)

  # window definition needs to be propagated into merge_windows() funtion
  starts <- genomic_bins_starts(start = 1, end = sequence_full_length, bin_size = window_size)
  if (length(starts) == 1) {
    ends <- sequence_full_length
  } else {
    ends <- c((starts[2 : length(starts)] - 1), sequence_full_length)
  }
  
  scores <- foreach (i = seq_along(starts), .combine = c, .export = c("log_messages", "extract_kmers", "seq_win_score_int")) %dopar% {
    seq_win_score_int(1, window_size, kmer, fasta_sequence[starts[i] : ends[i]], fraction_p)
  }
  if ((sum(is.na(scores))>0) || (length(scores) == 0)) {
    stop("sequence_window_score did not produce result")
  }
  if (log_messages != "") cat("\n", window_size, file = log_messages, append = TRUE)
  if (log_messages != "") cat("\n", length(fasta_sequence), file = log_messages, append = TRUE)

  return(scores)
}
