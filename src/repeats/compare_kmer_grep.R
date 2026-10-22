compare_kmer_grep = function(sequence_kmers, sequence_to_realign, max_size_dif, string_length, kmer = 6) {

  realign_string_length <- nchar(sequence_to_realign)
  if(!(realign_string_length %in% floor(string_length * (1 - max_size_dif)) : ceiling(string_length * (1 + max_size_dif)))) {
    return(sequence_to_realign)
  }
  sequence_to_realign_rev <- rev_comp_string(sequence_to_realign)
  score_orientation <- function(candidate) {
    candidate_extended <- paste0(rep(candidate, 2), collapse = "")
    candidate_kmers <- vapply(seq_len(realign_string_length), function(position) {
      substr(candidate_extended, position, position + kmer - 1)
    }, character(1))
    vapply(seq_len(realign_string_length), function(shift) {
      shifted_kmers <- candidate_kmers[((seq_len(length(sequence_kmers)) + shift - 2) %% realign_string_length) + 1]
      sum(sequence_kmers == shifted_kmers)
    }, numeric(1))
  }

  scores <- score_orientation(sequence_to_realign)
  scores_rev <- score_orientation(sequence_to_realign_rev)
  if (max(scores, scores_rev) == 0) return(sequence_to_realign)
  if (max(scores) > max(scores_rev)) {
    shifted <- paste0(strsplit(sequence_to_realign, split = "")[[1]][((seq_len(realign_string_length) + which.max(scores) - 2) %% realign_string_length) + 1], collapse = "")
  } else {
    shifted <- paste0(strsplit(sequence_to_realign_rev, split = "")[[1]][((seq_len(realign_string_length) + which.max(scores_rev) - 2) %% realign_string_length) + 1], collapse = "")
  }
  
  remove(sequence_kmers, sequence_to_realign, sequence_to_realign_rev)
  # gc()
  return(shifted)
}
# class_sequence = "aagcagtttcacagatagcttctttctagtttttatctggggatattcggtttttccccataggcctcaatgggctcccaaatgtcccttcgcagattctccaaaagagtgtttccaacctgctgaatcaaaagaaaggtttaactctgtgagatgaatccacacatcaca"
# sequence_to_realign = tolower("TCTAGTTTTTATCTGAAGATATTTTCTTTTTCACCATAGGCCTCAATGAGCTCCCAAATGTCCTTTCGCAGATTCTACAAAAACAGTGTTTCCAAACTGCTGAATCAAAAGAAAGGTTTAACTCTGTGAGATGAATGCACACATCACAAAGCAGTTTCTCAGAAAGCTTCCT")
