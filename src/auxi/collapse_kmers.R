collapse_kmers <- function(kmer_counts, kmer_names, max_edit = 3, verbose = FALSE) {

  kmer_names <- kmer_names[order(kmer_counts, decreasing = FALSE)]
  kmer_counts <- kmer_counts[order(kmer_counts, decreasing = FALSE)]
  return_list = list()
  if(max_edit == 0) {
    for(i in seq_along(kmer_counts)) {
      return_list[[i]] <- list(count = kmer_counts[i], kmers = kmer_names[[i]])
    }
    remove(kmer_counts, kmer_names, max_edit)
    gc()
    return(return_list)
  }

  collapsed_kmers_counts <- NULL
  collapsed_kmer_names <- list()
  i <- 1
  while (i <= length(kmer_counts)) {
    if(verbose) print(paste0(i, "/", length(kmer_counts)))
    collapsed_kmers_counts <- c(collapsed_kmers_counts, kmer_counts[i])
    collapsed_kmer_names <- append(collapsed_kmer_names, list(kmer_names[i]))

    if (i < length(kmer_counts)) {
      distances <- adist(kmer_names[i], kmer_names[(i+1):length(kmer_counts)], costs = list(insertions = 6, deletions = 6, substitutions = 1)) 
      merge_kmers <- NULL
      merge_kmers <- which(distances <= max_edit) + i
      collapsed_kmers_counts[length(collapsed_kmers_counts)] <- collapsed_kmers_counts[length(collapsed_kmers_counts)] + sum(kmer_counts[merge_kmers])
      collapsed_kmer_names[length(collapsed_kmer_names)][[1]] <- c(collapsed_kmer_names[length(collapsed_kmer_names)][[1]], kmer_names[merge_kmers])
      if(length(merge_kmers) != 0) {
        kmer_counts <- kmer_counts[-merge_kmers]
        kmer_names <- kmer_names[-merge_kmers]
      }
    }
    i = i + 1
  }
  
  for (i in seq_along(collapsed_kmers_counts)) {
    return_list[[i]] <- list(count = collapsed_kmers_counts[i], kmers = collapsed_kmer_names[[i]])
  }
  remove(collapsed_kmers_counts, collapsed_kmer_names, kmer_names, kmer_counts, max_edit, distances)
  gc()
  return(return_list)
}
