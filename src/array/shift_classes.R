shift_classes <- function(arrays, kmer = 6) {
  # This takes part of the array df with the same class to shift the representative to align with the one that's most common
  # start	end	seqID	numID	score	top_N	top_5_N	representative	class

  arrays$width <- arrays$end - arrays$start
  arrays$importance <- arrays$width * arrays$score
  if(sum(which(arrays$score == 0)) > 0) arrays$importance[arrays$score == 0] <- 1 * arrays$width[arrays$score == 0]

  class_sequence <- arrays$representative[which.max(arrays$importance)]
  string_length = nchar(class_sequence)
  class_sequence_extended = paste(rep(class_sequence, ceiling((string_length + kmer) / string_length)), collapse = "")
  class_sequence_kmers = unlist(lapply(seq_len(string_length), function(X) substr(class_sequence_extended, X, (X + kmer - 1))))

  sequences_to_realign <- arrays$representative
  sequences_shifted = vector(mode = "character", length = length(sequences_to_realign))
  for (i in seq_along(sequences_to_realign)) {
    sequences_shifted[i] = compare_kmer_grep(class_sequence_kmers, sequences_to_realign[i], 1, string_length, kmer = kmer)
  }
  remove(arrays, class_sequence, string_length, class_sequence_extended, class_sequence_kmers, sequences_to_realign)
  gc()
  return(sequences_shifted)
}