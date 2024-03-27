sequence_window_score <- function(fasta_sequence, window_size, kmer = 10) {
  # TODO check if changing these settings below can make the script work better,
  # although these were optimised
  fraction_p <- 0.5
  if ((window_size / 2) <= kmer) stop("sequence_window_score: window size is too small")

  sequence_full_length <- length(fasta_sequence)

  # window definition needs to be propagated into merge_windows() funtion
  sliding_window_distance <- ceiling(window_size / 4)
  if(sliding_window_distance >= sequence_full_length) {
    starts <- 1
    ends <- sequence_full_length
  } else {
    starts <- genomic_bins_starts(start = 1, end = sequence_full_length, bin_size = sliding_window_distance)
    starts <- starts[starts < sequence_full_length]
    if (length(starts) == 1) {
      ends <- sequence_full_length
    } else {
      ends <- c((starts[2 : length(starts)] - 1), sequence_full_length) + window_size # This makes overlapping windows!
    }
    ends[ends > sequence_full_length] <- sequence_full_length
  }
  
  scores <- foreach (i = seq_along(starts), .combine = c, .export = c("extract_kmers", "seq_win_score_int")) %dopar% {
    seq_win_score_int(1, window_size, kmer, fasta_sequence[starts[i] : ends[i]], fraction_p)
  }
  if ((sum(is.na(scores))>0) || (length(scores) == 0)) {
    stop("sequence_window_score did not produce valid result")
  }
  return(scores)
}
