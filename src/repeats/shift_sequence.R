shift_sequence <- function(sequence, kmer = 4) {
  string_length <- nchar(sequence)
  sequence <- paste0(rep(sequence, ceiling((string_length + kmer) / string_length)), collapse = "")
  kmers <- unlist(lapply(seq_len(string_length), function(X) substr(sequence, X, (X + kmer - 1))))
  kmers_scores <- unlist(lapply(kmers, kmer_hash_score))

  kmers <- c(kmers, kmers)

  hash_values <- unlist(lapply(X = seq_len(string_length), hash_sequence(kmers[X : (X + string_length - 1)])))

  return(sequence)
}