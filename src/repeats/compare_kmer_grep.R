compare_kmer_grep = function(sequence_kmers, sequence_to_realign, max_size_dif, string_length, kmer = 6) {

  realign_string_length <- nchar(sequence_to_realign)
  if(!(realign_string_length %in% floor(string_length * (1 - max_size_dif)) : ceiling(string_length * (1 + max_size_dif)))) {
    return(sequence_to_realign)
  }
  seq_to_realign_kmers_ext = paste0(rep(sequence_to_realign, 1 + ceiling((realign_string_length + kmer) / realign_string_length)), collapse = "")
  seq_to_realign_kmers = unlist(lapply(seq_len(realign_string_length * ceiling((realign_string_length + kmer) / realign_string_length)), function(X) substr(seq_to_realign_kmers_ext, X, (X + kmer - 1))))

  kmers_distances = rep(0, length(sequence_kmers))
  for(i in seq_along(sequence_kmers)) {
    min = grep(sequence_kmers[i], seq_to_realign_kmers[i : length(seq_to_realign_kmers)])
    if(length(min) == 0) min = 0 else min = min(min)
    kmers_distances[i] = min
  }
  if(sum(kmers_distances) == 0) return(sequence_to_realign)
  shift = as.numeric(names(table(kmers_distances[kmers_distances != 0]))[which.max(table(kmers_distances[kmers_distances != 0]))])
  shifted = paste(unlist(lapply(seq_to_realign_kmers[shift : (shift + realign_string_length - 1)], function(X) strsplit(X, split = "")[[1]][1])), collapse = "")
  remove(sequence_kmers, sequence_to_realign, seq_to_realign_kmers_ext, seq_to_realign_kmers, shift, kmers_distances)
  gc()
  return(shifted)
}
# class_sequence = "aagcagtttcacagatagcttctttctagtttttatctggggatattcggtttttccccataggcctcaatgggctcccaaatgtcccttcgcagattctccaaaagagtgtttccaacctgctgaatcaaaagaaaggtttaactctgtgagatgaatccacacatcaca"
# sequence_to_realign = tolower("TCTAGTTTTTATCTGAAGATATTTTCTTTTTCACCATAGGCCTCAATGAGCTCCCAAATGTCCTTTCGCAGATTCTACAAAAACAGTGTTTCCAAACTGCTGAATCAAAAGAAAGGTTTAACTCTGTGAGATGAATGCACACATCACAAAGCAGTTTCTCAGAAAGCTTCCT")
