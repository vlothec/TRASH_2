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

  scores <- vector(mode = "numeric", length = length(starts))

  for (i in seq_along(starts)) {
    kmers = unlist(lapply(X = starts[i] : (ends[i] - kmer), FUN = extract_kmers, kmer, fasta_sequence))
    counts_kmers <- table(kmers)
    counts_kmers <- counts_kmers[!grepl("n", names(counts_kmers))]
    counts_kmers <- counts_kmers[!grepl("N", names(counts_kmers))]

    if (length(counts_kmers) < kmer) {
      scores[i] <- 100
      next
    }
    total_kmers <- sum(counts_kmers)
    counts_kmers <- sort(counts_kmers, decreasing = TRUE)

    counts_kmers = unlist(lapply(seq_along(counts_kmers), FUN = function(x, vector) return(sum(vector[1 : x])), counts_kmers))
    scores[i] <- 100 * (min(which(counts_kmers / total_kmers > fraction_p)) - 1) / (total_kmers) / fraction_p
  }
  return(scores)
}