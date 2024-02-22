kmer_hash_score = function(kmer) {
  kmer = strsplit(kmer, split = "")[[1]]
  actg_val = 0 * (hash_single[j] == "A") + 1 * (hash_single[j] == "C") + 2 * (hash_single[j] == "T") + 3 * (hash_single[j] == "G") 
  hash_val[i] = hash_val[i] + (4 ^ j) * actg_val
  return(0)
}