seq_win_score_int <- function(start, end, kmer, fasta_extraction, fraction_p) {
  kmers <- unlist(lapply(X = start : (end - kmer), FUN = extract_kmers, kmer, fasta_extraction))
  counts_kmers <- table(kmers)
  counts_kmers <- counts_kmers[!grepl("n", names(counts_kmers))]
  counts_kmers <- counts_kmers[!grepl("N", names(counts_kmers))]

  total_kmers <- sum(counts_kmers)
  if (total_kmers < (kmer * 2)) return(100)

  counts_kmers <- sort(counts_kmers, decreasing = TRUE)
  counts_kmers <- unlist(lapply(seq_along(counts_kmers), FUN = function(x, vector) return(sum(vector[1 : x])), counts_kmers))
  score <- 100 * (min(which(counts_kmers / total_kmers > fraction_p)) - 1) / (total_kmers) / fraction_p
  remove(start, end, kmer, fasta_extraction, fraction_p, counts_kmers, total_kmers)
  return(score)
}
