split_and_check_arrays <- function(start, end, sequence, seqID, numID) {
  
  kmer <- 12
  min_kmers_count = 2
  
  window_size = 1000 #small windows which are compared in a row to find array break points
  max_window_count_percentage_to_score = 0.2
  score_weight_similar_count = 1
  score_weight_similar_top_distance = 1
  
  array_break_max_threshold = 1 # out of 6
  
  small_window_for_N_count = 1000 
  small_window_step_for_N_count = 100
  small_window_min_percentage_of_distances = small_window_for_N_count / 10
  

  kmers_list <- unlist(lapply(X = (start : (end - kmer)), FUN = extract_kmers, kmer, sequence))
  counts_kmers <- table(kmers_list)
  kmer_names <- names(counts_kmers)

  #clean up N containing kmers
  counts_kmers <- counts_kmers[!grepl("N", kmer_names)]
  kmer_names <- kmer_names[!grepl("N", kmer_names)]
  counts_kmers <- counts_kmers[!grepl("n", kmer_names)]
  kmer_names <- kmer_names[!grepl("n", kmer_names)]

  if (length(counts_kmers) == 0) {
    return(NULL)
  }

  kmer_names = kmer_names[counts_kmers >= min_kmers_count]
  counts_kmers = counts_kmers[counts_kmers >= min_kmers_count]

  if (length(counts_kmers) == 0) {
    return(NULL)
  }

  collapsed_kmers <- collapse_kmers(counts_kmers, kmer_names)

  for (i in seq_along(collapsed_kmers)) {
    collapsed_kmers[[i]]$locations <- which(kmers_list %in% collapsed_kmers[[i]]$kmers)
  }

  plot(x = NULL, y = NULL, xlim = c(0,100000), ylim = c(0, length(collapsed_kmers)))
  for (i in seq_along(collapsed_kmers)) {
   points(x = collapsed_kmers[[i]]$locations, y = rep(i, collapsed_kmers[[i]]$count), pch = ".")
  }
  
  window_starts = NULL
  window_ends = NULL
  
  window_starts <- genomic_bins_starts(start = start, end = end, bin_size = window_size)
  if (length(window_starts) < 2) {
    window_ends <- end
  } else {
    window_ends <- c((window_starts[2:length(window_starts)] - 1), end)
  }
  if (length(window_ends) != length(window_starts)) window_ends <- end
  
  window_top_1_distance <- rep(0, length(window_starts))
  window_top_2_distance <- rep(0, length(window_starts))
  window_top_3_distance <- rep(0, length(window_starts))
  window_top_4_distance <- rep(0, length(window_starts))
  window_top_5_distance <- rep(0, length(window_starts))
  window_count <- rep(0, length(window_starts))
  
  distances <- NULL
  kmer_starts <- NULL
  for (i in seq_along(collapsed_kmers)) {
    collapsed_kmers[[i]]$distances = (collapsed_kmers[[i]]$locations[2:length(collapsed_kmers[[i]]$locations)] -
                                        collapsed_kmers[[i]]$locations[1:(length(collapsed_kmers[[i]]$locations)-1)])
    
    distances <- c(distances, collapsed_kmers[[i]]$distances)
    kmer_starts <- c(kmer_starts, collapsed_kmers[[i]]$locations[1:(length(collapsed_kmers[[i]]$locations)-1)])
  }
  
  for (i in seq_along(window_starts)) {
    window_count[i] <- length(which(kmer_starts >= window_starts[i] & kmer_starts <= window_ends[i]))
    window_distances <- table(distances[which(kmer_starts >= window_starts[i] & kmer_starts <= window_ends[i])])
    window_distances <- window_distances[order(window_distances, decreasing = TRUE)]

    window_top_1_distance[i] <- if(length(as.numeric(names(window_distances)[1])) == 0) 0 else as.numeric(names(window_distances)[1]) 
    window_top_2_distance[i] <- if(length(as.numeric(names(window_distances)[2])) == 0) 0 else as.numeric(names(window_distances)[2]) 
    window_top_3_distance[i] <- if(length(as.numeric(names(window_distances)[3])) == 0) 0 else as.numeric(names(window_distances)[3]) 
    window_top_4_distance[i] <- if(length(as.numeric(names(window_distances)[4])) == 0) 0 else as.numeric(names(window_distances)[4]) 
    window_top_5_distance[i] <- if(length(as.numeric(names(window_distances)[5])) == 0) 0 else as.numeric(names(window_distances)[5]) 
  }
  
  connected_next_windows <- rep(0, length(window_starts))
  
  for(i in 1 : (length(window_starts) - 1)) {
    next_distances = c(window_top_1_distance[i + 1], window_top_2_distance[i + 1], window_top_3_distance[i + 1], window_top_4_distance[i + 1], window_top_5_distance[i + 1])
    next_distances = next_distances[!is.na(next_distances)]
    next_distances = next_distances[next_distances != 0]
    if(length(next_distances) == 0) next
    connected_next_windows[i] <- connected_next_windows[i] + score_weight_similar_count * ((abs(window_count[i] - window_count[i + 1]) / (max(window_count[i], window_count[i + 1]) + 1)) < max_window_count_percentage_to_score)
    connected_next_windows[i] <- connected_next_windows[i] + score_weight_similar_top_distance * (window_top_1_distance[i] %in% next_distances)
    connected_next_windows[i] <- connected_next_windows[i] + score_weight_similar_top_distance * (window_top_2_distance[i] %in% next_distances)
    connected_next_windows[i] <- connected_next_windows[i] + score_weight_similar_top_distance * (window_top_3_distance[i] %in% next_distances)
    connected_next_windows[i] <- connected_next_windows[i] + score_weight_similar_top_distance * (window_top_4_distance[i] %in% next_distances)
    connected_next_windows[i] <- connected_next_windows[i] + score_weight_similar_top_distance * (window_top_5_distance[i] %in% next_distances)
  }
  connected_next_windows[length(connected_next_windows)] = score_weight_similar_count + 5 * score_weight_similar_top_distance
  
  array_breaks_estimate_starts = window_starts[connected_next_windows <= array_break_max_threshold]
  array_breaks_estimate_ends = window_ends[connected_next_windows <= array_break_max_threshold]
  
  array_breaks_estimate_starts_temp = array_breaks_estimate_starts
  array_breaks_estimate_starts = array_breaks_estimate_starts[!((array_breaks_estimate_starts-1) %in% array_breaks_estimate_ends)]
  array_breaks_estimate_ends = array_breaks_estimate_ends[!((array_breaks_estimate_ends+1) %in% array_breaks_estimate_starts_temp)]

  
  window_starts <- genomic_bins_starts(start = start, end = end, bin_size = small_window_step_for_N_count)
  if (length(window_starts) < 2) {
    window_ends <- end
  } else {
    window_ends <- c((window_starts[2:length(window_starts)] - 1), end)
  }
  if (length(window_ends) != length(window_starts)) window_ends <- end
  window_ends = window_ends - small_window_step_for_N_count + small_window_for_N_count
  window_ends[window_ends > end] = end
  
  moving_top_distance = vector(length = length(window_starts), mode = "numeric")
  for(i in seq_along(window_starts)) {
    which_distances = (kmer_starts >= window_starts[i] & kmer_starts <= window_ends[i])
    if(sum(which_distances) > small_window_min_percentage_of_distances) {
      moving_top_distance[i] = as.numeric(names(which.max(table(distances[which_distances]))))
    }
  }
  
  plot(x = window_starts[moving_top_distance != 0], y = moving_top_distance[moving_top_distance != 0], pch = 16, ylim = c(0, max(moving_top_distance)))
  abline(v = c(34450,39110,44870))
  table(moving_top_distance[moving_top_distance != 0])

  collapsed_kmers_counts = NULL
  collapsed_kmers_most_common_distance = NULL
  for(i in seq_along(collapsed_kmers)) {
    collapsed_kmers_counts = c(collapsed_kmers_counts, collapsed_kmers[[i]]$count)
    collapsed_kmers[[i]]$most_common_distance = as.numeric(names(which.max(table(collapsed_kmers[[i]]$distances))))
    collapsed_kmers_most_common_distance = c(collapsed_kmers_most_common_distance, collapsed_kmers[[i]]$most_common_distance )
  }
  
  
  
  return(arrays)
}
#Test
#test_fasta = read_fasta_and_list("../../testing_fastas/SUPER_15_extraction.fasta")
#start = 1
#end = 100000
#sequence = test_fasta[[1]]
#window_size = 1000
#split_and_check_arrays(1,10000,test_fasta[[1]], "CP068268", 1)

# Split arrays:kmer distances, use colapsed kmers;
# start with most common kmer, estimate its range
# (consider locations that are reatively nearby);
# do the same with each kmer; overlap the kmer ranges;
# find a bg_range that has the most overlapping ranges;
# subtract ranges that overlap with the big_range;
# identify next big_range and so on until no more ranges left;
# using kmers(BIG_range_1) %in% kmers(BIG_range_1) and
# by comparing most frequent N, check if they contain
# the same repeat; if all BIG_ranges contain the same
# repeat, keep as **simple region**; if there are different
# repeats, use the initial ranges to split the region into
# two (or more), but give each of them an overlap (at least
# N*10 but not more than 10% of each region) so that the
# edge repeats can be properly identified, assign them as
# **simple region**; if there are BIG_ranges with more than
# one repeat, but each of them appearing multiple times (like
# BIG_range_1, BIG_range_2, BIG_range_1, BIG_range_2), then
# assign as a **complex region**

# colapsed kmers:
# start at the most common kmer ->
# identify all other kmers that are
# in edit distance 1 to that kmer,
# merge their locations; move to the
# next un-merged kmer until the list is over