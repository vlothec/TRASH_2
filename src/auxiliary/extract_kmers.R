extract_kmers <- function(X, kmer, fasta_sequence) {
  end = X + kmer - 1
  end = length(fasta_sequence) * (end > length(fasta_sequence)) + end * (end <= length(fasta_sequence))
  return(paste(fasta_sequence[X : end], collapse = "")) # nolint
}
