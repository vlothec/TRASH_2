extract_kmers <- function(X, kmer, fasta_sequence) { # I've replaced this function in the script and is liekly not executed anywhere, passing the huge sequence into in multiple times seems problematic
  end = X + kmer - 1
  return(paste(fasta_sequence[X : end], collapse = "")) # nolint
}
