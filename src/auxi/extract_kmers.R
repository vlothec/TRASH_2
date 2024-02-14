extract_kmers <- function(X, kmer, fasta_sequence) {
  end = X + kmer - 1
  return(paste(fasta_sequence[X : end], collapse = "")) # nolint
}
