split_and_check_arrays <- function(start, end, sequence, seqID, numID, max_repeat) {
  
  ### Extract kmers and calculate distances ===============================================================
  kmer <- 12
  min_kmers_count = 2
  
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
  
  distances <- NULL
  kmer_starts <- NULL
  for (i in seq_along(collapsed_kmers)) {
    collapsed_kmers[[i]]$distances = (collapsed_kmers[[i]]$locations[2:length(collapsed_kmers[[i]]$locations)] -
                                        collapsed_kmers[[i]]$locations[1:(length(collapsed_kmers[[i]]$locations)-1)])
    kmer_starts <- c(kmer_starts, collapsed_kmers[[i]]$locations[1:(length(collapsed_kmers[[i]]$locations)-1)])
    distances <- c(distances, collapsed_kmers[[i]]$distances)
    collapsed_kmers[[i]]$distances = collapsed_kmers[[i]]$distances[collapsed_kmers[[i]]$distances <= max_repeat]
  }
  ### =====================================================================================================


  ### Find breaks  ========================================================================================
  window_size = max_repeat * 2
  if(window_size < 1000) window_size = 1000
  window_step = 50
  window_starts = NULL
  window_ends = NULL
  window_ends_compare = NULL
  min_windows_comparison_score_to_detach_array = 0.2
  min_windows_comparison_score_to_split_array = 0.1
  windows_comparison_score = NULL
  array_overlaps = ceiling(max_repeat/2)
  array_breaks_coordinates = NULL
  
  if((end - window_size) > start) { # if it's not, windows_comparison_score is not set and arrays are not split
    window_starts <- genomic_bins_starts(start = start, end = (end - window_size), bin_size = window_step)
    if (length(window_starts) < 2) {
      window_ends <- (end - window_size)
    } else {
      window_ends <- c((window_starts[2:length(window_starts)] - 1), (end - window_size))
    }
    if (length(window_ends) != length(window_starts)) window_ends <- (end - window_size)
    window_ends = window_ends + window_size - window_step
    window_ends[window_ends > (end - window_size)] = (end - window_size)
    
    window_starts_compare = window_ends + 1
    
    if (length(window_starts_compare) < 2) {
      window_ends_compare <- end
    } else {
      window_ends_compare <- c((window_starts_compare[2:length(window_starts_compare)] - 1), end)
    }
    if (length(window_ends_compare) != length(window_starts_compare)) window_ends_compare <- end
    window_ends_compare = window_ends_compare + window_size - window_step
    window_ends_compare[window_ends_compare > end] = end
    
    windows_comparison_score <- rep(0, length(window_starts))
    
    for(i in seq_along(window_starts_compare)) {
      kmers_window_A = kmers_list[window_starts[i] : window_ends[i]]
      kmers_window_B = kmers_list[window_starts_compare[i] : window_ends_compare[i]]
      windows_comparison_score[i] = sum(kmers_window_A %in% kmers_window_B) / window_size
    }
  }
  
  if (length(windows_comparison_score) == 0) {
    # the region is too small to compare anything, analyse it as is
    arrays = data.frame(start = start, end = end, seqID = seqID, numID = numID, windows_comparison_score = -1)
  } else if(sum(windows_comparison_score > min_windows_comparison_score_to_detach_array) == 0) {
    # there is no windows above the threshold to split arrays, likely a mislabeled region or a complex one or contains
    # repeats longer than the settings allow for
    arrays = data.frame(start = start, end = end, seqID = seqID, numID = numID, score = mean(windows_comparison_score))
  } else {
    #attempt splitting and save each array, ALLOW overlaps equal to 1 max repeat size 
    window_event = "new_top"
    window_event_position = 1
    i = which(windows_comparison_score > min_windows_comparison_score_to_detach_array)[1] + 1
    while (i <= length(window_starts_compare)) {
      while (i <= length(window_starts_compare) & window_event[length(window_event)] == "new_top") {
        if(windows_comparison_score[i] <= min_windows_comparison_score_to_split_array) {
          window_event = c(window_event, "new_bottom")
          window_event_position = c(window_event_position , i)
        }
        i = i + 1
      }
      while (i <= length(window_starts_compare) & window_event[length(window_event)] == "new_bottom") {
        if(windows_comparison_score[i] >= min_windows_comparison_score_to_detach_array) {
          window_event = c(window_event, "new_top")
          window_event_position = c(window_event_position , i)
        }
        i = i + 1
      }
    }
    if(length(window_event) == 1) {
      # no second array found, report as is
      arrays = data.frame(start = start, end = end, seqID = seqID, numID = numID, score = mean(windows_comparison_score))
    } else {
      # there are multiple arrays, find them all
      array_breaks = NULL
      for (i in 1 : (length(window_event) - 1)) {
        if(window_event[i] == "new_bottom" & window_event[i+1] == "new_top") {
          array_breaks = c(array_breaks, (window_event_position[i] + which.min(windows_comparison_score[window_event_position[i] : window_event_position[i+1]]) - 1))
        }
      }
      array_breaks_coordinates = window_ends[array_breaks] + 1
      arrays = data.frame(start = 1, 
                          end = (array_breaks_coordinates[1] - 1 + array_overlaps), 
                          seqID = seqID, 
                          numID = numID, 
                          score = mean(windows_comparison_score[1 : array_breaks[1]])) # This score is not perfect (especially with short arrays), don't use too much
      if(length(array_breaks) > 1) {
        for(i in 1 : (length(array_breaks) - 1)) {
          arrays = rbind(arrays, list(start = array_breaks_coordinates[i] - array_overlaps,
                                      end = (array_breaks_coordinates[i+1] - 1 + array_overlaps),
                                      seqID = seqID, 
                                      numID = numID, 
                                      score = mean(windows_comparison_score[array_breaks[i] : array_breaks[i + 1]])))
        }
      }
      arrays = rbind(arrays, list(start = array_breaks_coordinates[length(array_breaks)],
                                  end = end,
                                  seqID = seqID, 
                                  numID = numID, 
                                  score = mean(windows_comparison_score[array_breaks[length(array_breaks)] : length(windows_comparison_score)])))
    }
  }
  if(!inherits(arrays, "data.frame")) { # sanity check, this should not happen
    arrays = data.frame(start = start, end = end, seqID = seqID, numID = numID, score = -1)
    print("Check array identification")
  } else if (nrow(arrays) == 0) {
    arrays = data.frame(start = start, end = end, seqID = seqID, numID = numID, score = -1)
    print("Check array identification")
  }
  ### =====================================================================================================
  
  
  ### For each array  =====================================================================================
  arrays$top_N = 0
  arrays$top_5_N = "."
  arrays$representative = ""
  for(i in seq_len(nrow(arrays))) {
    ## recalculate distances if more than 1 array =============================================
    if(nrow(arrays) != 1) {
      
      kmers_list <- unlist(lapply(X = (arrays$start[i] : (arrays$end[i] - kmer)), FUN = extract_kmers, kmer, sequence))
      counts_kmers <- table(kmers_list)
      kmer_names <- names(counts_kmers)
      #clean up N containing kmers
      counts_kmers <- counts_kmers[!grepl("N", kmer_names)]
      kmer_names <- kmer_names[!grepl("N", kmer_names)]
      counts_kmers <- counts_kmers[!grepl("n", kmer_names)]
      kmer_names <- kmer_names[!grepl("n", kmer_names)]
      
      if (length(counts_kmers) == 0) {
        print("Empty array, should not happen")
        next
      }
      
      kmer_names = kmer_names[counts_kmers >= min_kmers_count]
      counts_kmers = counts_kmers[counts_kmers >= min_kmers_count]
      if (length(counts_kmers) == 0) {
        print("Empty array, should not happen")
        next
      }
      
      collapsed_kmers <- collapse_kmers(counts_kmers, kmer_names)
      for (j in seq_along(collapsed_kmers)) {
        collapsed_kmers[[j]]$locations <- which(kmers_list %in% collapsed_kmers[[j]]$kmers)
      }
      
      distances <- NULL
      kmer_starts <- NULL
      for (j in seq_along(collapsed_kmers)) {
        collapsed_kmers[[j]]$distances = (collapsed_kmers[[j]]$locations[2:length(collapsed_kmers[[j]]$locations)] -
                                            collapsed_kmers[[j]]$locations[1:(length(collapsed_kmers[[j]]$locations)-1)])
        kmer_starts <- c(kmer_starts, collapsed_kmers[[j]]$locations[1:(length(collapsed_kmers[[j]]$locations)-1)])
        distances <- c(distances, collapsed_kmers[[j]]$distances)
        collapsed_kmers[[j]]$distances = collapsed_kmers[[j]]$distances[collapsed_kmers[[j]]$distances <= max_repeat]
      }
      kmer_starts = kmer_starts - 1 + arrays$start[i]
    }
    ## ========================================================================================
    
    ## Find N using kmer distances ============================================================
    small_window_for_N_count = 1000 
    small_window_step_for_N_count = 100
    small_window_min_percentage_of_distances = small_window_for_N_count / 10
    window_starts <- genomic_bins_starts(start = arrays$start[i], end = arrays$end[i], bin_size = small_window_step_for_N_count)
    
    if (length(window_starts) < 2) {
      window_ends <- arrays$end[i]
    } else {
      window_ends <- c((window_starts[2:length(window_starts)] - 1), arrays$end[i])
    }
    if (length(window_ends) != length(window_starts)) window_ends <- arrays$end[i]
    window_ends = window_ends - small_window_step_for_N_count + small_window_for_N_count
    window_ends[window_ends > arrays$end[i]] = arrays$end[i]
    
    moving_top_distance = vector(length = length(window_starts), mode = "numeric")
    
    for(j in seq_along(window_starts)) {
      which_distances = (kmer_starts >= window_starts[j] & kmer_starts <= window_ends[j])
      if(sum(which_distances) > small_window_min_percentage_of_distances) {
        moving_top_distance[j] = as.numeric(names(which.max(table(distances[which_distances]))))
      }
    }
    top_N_array = sort(table(moving_top_distance[moving_top_distance != 0]), decreasing = TRUE)

    # merge N values that are only 1 bp apart, 
    # TODO: consider removing low count N values in case a region contains them all, 
    # or limit to top X Ns to be merged only
    # or limit to Ns with at least a min count to be allowed to be merged
    top_N_array = top_N_array[order(as.numeric(names(top_N_array)))]
    j = 1
    last_name = as.numeric(names(top_N_array)[j])
    while(j < length(top_N_array)) {
      if(as.numeric(names(top_N_array)[j + 1]) == (last_name + 1)) {
        last_name = as.numeric(names(top_N_array)[j + 1])
        top_N_array[j] = top_N_array[j] + top_N_array[j + 1]
        names(top_N_array)[j] = paste0(names(top_N_array)[j], ",", as.character(last_name))
        top_N_array = top_N_array[-(j+1)]
      } else {
        last_name = as.numeric(names(top_N_array)[j + 1])
        j = j + 1
      }
    }
    top_N_array = sort(top_N_array, decreasing = TRUE)
    
    top_N_distances = NULL
    count_Ns = length(top_N_array)
    if(count_Ns != 0) {
      if(count_Ns >= 5) count_Ns = 5
      for (j in seq_len(count_Ns)) {
        arrays$top_5_N[i] = paste0(arrays$top_5_N[i], "N_", names(top_N_array[j]), "_Count_", top_N_array[j], ".")
      }
      top_N_distances = strsplit(arrays$top_5_N[i], split = "[.]")[[1]][2]
      top_N_distances = strsplit(top_N_distances, split = "[_]")[[1]][2]
      top_N_distances = strsplit(top_N_distances, split = "[,]")[[1]]
      top_N_distances = as.numeric(top_N_distances)
      arrays$top_N[i] = floor(mean(top_N_distances))
    }
    ## ========================================================================================
    
    ## Identify kmers likely forming the repeat ===============================================
    # use kmer_starts and distances, find the best max_repeat size window
    # then, within that window, find the best arrays$top_N[i] size window
    kmer_hist = hist(kmer_starts[distances %in% top_N_distances], 
                     breaks = genomic_bins_starts(arrays$start[i], (arrays$end[i] + 2 * max_repeat), bin_size = max_repeat),
                     plot = FALSE)
    window_with_representative_start = kmer_hist$breaks[which.max(kmer_hist$counts)]
    window_with_representative_end = window_with_representative_start + max_repeat - 1
    
    window_with_rep_kmer_starts = kmer_starts[kmer_starts >= window_with_representative_start & kmer_starts <= window_with_representative_end]
    window_with_rep_distances = distances[kmer_starts >= window_with_representative_start & kmer_starts <= window_with_representative_end]
    
    window_with_rep_kmer_starts = window_with_rep_kmer_starts[window_with_rep_distances %in% top_N_distances]
    window_with_rep_distances = window_with_rep_distances[window_with_rep_distances %in% top_N_distances]
    
    window_with_rep_kmer_starts = sort(window_with_rep_kmer_starts, decreasing = FALSE)
    
    window_hist = hist(window_with_rep_kmer_starts, 
                       breaks = genomic_bins_starts(window_with_representative_start, (window_with_representative_end + arrays$top_N[i]), bin_size = arrays$top_N[i]),
                       plot = FALSE)
    
    representative_start = window_hist$breaks[which.max(window_hist$counts)]
    representative_end = representative_start + arrays$top_N[i] - 1
    
    arrays$representative[i] = paste(sequence[representative_start : representative_end], collapse = "")
    
    ## Identify kmers likely forming the repeat ===============================================
    # Use the best kmer and extract up to 100 repeats, align and get consensus
    collapsed_kmers_counts = NULL
    collapsed_kmers_most_common_distance = NULL
    for(i in seq_along(collapsed_kmers)) {
      collapsed_kmers_counts = c(collapsed_kmers_counts, collapsed_kmers[[i]]$count)
      collapsed_kmers[[i]]$most_common_distance = as.numeric(names(which.max(table(collapsed_kmers[[i]]$distances))))
      collapsed_kmers_most_common_distance = c(collapsed_kmers_most_common_distance, collapsed_kmers[[i]]$most_common_distance)
    }
    
    kmer_starts[distances %in% top_N_distances]
    distances[distances %in% top_N_distances]
    
    
    ## ========================================================================================
    
    ## Plot arrays ============================================================================
    if(T) {
      gc_calc_window = 1000
      
      par(mfrow = c(3,1), mar = c(4,4,1,1))
      # plot 1: moving top distance (N)
      plot(x = window_starts[moving_top_distance != 0], y = moving_top_distance[moving_top_distance != 0], pch = 16, 
           ylim = c(0, max_repeat),
           xlim = c(start, end))
      if(length(array_breaks_coordinates) != 0) abline(v = array_breaks_coordinates)
      # plot 2: windows comparison score
      plot(x = window_starts_compare[window_starts_compare > arrays$start[i] & window_starts_compare < arrays$end[i]], 
           y = windows_comparison_score[window_starts_compare > arrays$start[i] & window_starts_compare < arrays$end[i]], pch = 16, 
           ylim = c(0, 1),
           xlim = c(start, end))
      if(length(array_breaks_coordinates) != 0) abline(v = array_breaks_coordinates)
      # plot 3: GC
      window_starts <- genomic_bins_starts(start = arrays$start[i], end = arrays$end[i], bin_size = gc_calc_window)
      gc_track = calculate.GC.in.windows(windows.starts = (window_starts - window_starts[1] + 1), sequence = sequence[arrays$start[i] : arrays$end[i]], bin.size = gc_calc_window)
      plot(x = window_starts, y = gc_track, ylim = c(0,1), pch = 16, xlim = c(start, end))
      if(length(array_breaks_coordinates) != 0) abline(v = array_breaks_coordinates)
      par(mfrow = c(1,1))
    }
    ## ========================================================================================
  }
  ### =====================================================================================================
  return(arrays)
}
#Test
#test_fasta = read_fasta_and_list("../../testing_fastas/SUPER_15_extraction.fasta")
#test_fasta = read_fasta_and_list("../../testing_fastas/ath_Chr1_extraction.fasta")
#test_fasta = read_fasta_and_list("../../testing_fastas/Rbrevi_chr1_extr1_edited.fasta")
#start = 1
#sequence = test_fasta[[1]][45000:92000]
#end = length(sequence)
#window_size = 1000
#max_repeat = 800
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