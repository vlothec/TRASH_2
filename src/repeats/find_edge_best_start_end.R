find_edge_best_start_end <- function(sequence_potential, representative) { # two sequences are character vectors
  kmer = 12
  if (length(sequence_potential) <= kmer) return(0) #default case for right edge, returning 1 will discard new repeat in the higher function
  if (length(representative) <= kmer) return(0)
  new_end <- length(sequence_potential)

  representative_kmers <- unlist(lapply(1:(length(representative) - kmer), function(X) paste(representative[X : (X + kmer - 1)], collapse = "")))
  potential_kmers <- unlist(lapply(1:(length(sequence_potential) - kmer), function(X) paste(sequence_potential[X : (X + kmer - 1)], collapse = "")))
  found_kmers <- which(potential_kmers %in% representative_kmers)
  scores <- rep(0, length(sequence_potential))
  for (i in seq_along(scores)) {
    scores[i] <- sum(found_kmers <= i) 
  }
  division_score <- rep(0, length(sequence_potential))
  for (i in seq_along(division_score)) {
    division_score[i] <- ((length(unique(scores[1 : i]))) / i) / ((length(unique(scores[i : length(scores)]))) / (length(scores) - i + 1))
  }
  division_score <- division_score * (seq(0.25,1,length.out=length(division_score)))
  new_end <- which.max(division_score) + kmer
  if(sum(potential_kmers[1 : new_end] %in% representative_kmers) / new_end > 0.15) {
    return(new_end)
  } else {
    return(0)
  }

}


# representative <- rev(array_representative_rev)
# paste(sequence_potential, collapse = "")
# paste(representative, collapse = "")
