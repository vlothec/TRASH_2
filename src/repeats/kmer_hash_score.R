kmer_hash_score = function(kmer) {
  kmer = strsplit(kmer, split = "")[[1]]
  hash_val = 0
  for(i in seq_along(kmer)) {
    actg_val = 0 * (kmer[i] == "a") + 1 * (kmer[i] == "c") + 2 * (kmer[i] == "t") + 3 * (kmer[i] == "g") 
    hash_val = hash_val + (4 ^ i) * actg_val
  }
  return(hash_val)
}