sequence_window_score <- function(fasta_sequence, window_size) {
  # TODO make it in parallel, now it uses lapply
  # TODO check if changing these settings below can make the script work better,
  # although these were optimised
  kmer <- 12
  fraction_p <- 0.5

  sequence_full_length <- length(fasta_sequence)

  starts <- genomic_bins_starts(start = 1, end = sequence_full_length, bin_size = window_size)
  if (length(starts) < 2) {
    ends <- sequence_full_length
  } else {
    ends <- c((starts[2:length(starts)] - 1), sequence_full_length)
  }
  if (length(ends) != length(starts)) ends <- sequence_full_length

  scores <- foreach (i = seq_along(starts), .combine = c, .export = c("extract_kmers", "seq_win_score_int")) %dopar% {
    seq_win_score_int(1, window_size, kmer, fasta_sequence[starts[i] : ends[i]], fraction_p)
  }
  return(scores)
}