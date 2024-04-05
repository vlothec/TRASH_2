seq_win_score_int <- function(start, end, kmer, fasta_extraction, fraction_p) {
  if((end - start) <= kmer) return(100)
  # kmers <- unlist(lapply(X = start : (end - kmer), FUN = extract_kmers, kmer, fasta_extraction))
  kmers <- unlist(lapply(X = (start : (end - kmer)), function(X) return(paste(fasta_extraction[X : (X + kmer - 1)], collapse = ""))))
  counts_kmers <- table(kmers)
  counts_kmers <- counts_kmers[!grepl("n", names(counts_kmers))]
  counts_kmers <- counts_kmers[!grepl("N", names(counts_kmers))]

  if (sum(counts_kmers) < (kmer * 2)) return(100)
  score <- (100 * sum(counts_kmers[counts_kmers == 1]) / sum(counts_kmers))
  remove(start, end, kmer, fasta_extraction, fraction_p, counts_kmers, kmers)
  gc()
  return(score)

  # counts_kmers <- sort(counts_kmers, decreasing = TRUE)
  # counts_kmers <- unlist(lapply(seq_along(counts_kmers), FUN = function(x, vector) return(sum(vector[1 : x])), counts_kmers))
  # score <- 100 * (min(which(counts_kmers / total_kmers > fraction_p)) - 1) / (total_kmers) / fraction_p
  # remove(start, end, kmer, fasta_extraction, fraction_p, counts_kmers, total_kmers)
  # return(score)
}
# fasta_extraction = strsplit("actgctagtcgattcgaactgctaccgtcgatcgaactgctcttagtcgatcgaactctgctagtcgatcgattcccgactgctagtcgattcgaactgctaccgtcgatcgaactgagtcgatcgattcccgactgctagtcgattcgaactgctaccgtcgatcgaactgctcttagtcgatcgaactctgctagtcgatcgattcccg", split = "")[[1]]
# start = 1
# end = length(fasta_extraction)
# fraction_p = 0.5
# kmer = 10
