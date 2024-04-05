split_and_check_arrays <- function(start, end, sequence, seqID, numID, arrID, max_repeat, min_repeat, mafft, temp_dir, src_dir, sink_output, kmer = 10) {
  # sink(file.path(temp_dir, paste0(seqID, "_", arrID, "_logfile.txt")))

  ### Settings ===========================================================================================
  ## Find breaks
  window_step <- ceiling(max_repeat / 20)
  min_windows_comparison_score_to_detach_array <- 0.08
  min_windows_comparison_score_to_split_array <- 0.04
  array_overlaps <- 0 # setting just in case; overlap might be usefull for making sure edge repeats are mapped correctly
  # Filter kmers
  global_min_kmers_count <- 2 # obligatory filtering to decrease the number of unique kmers, removes single-copy ones
  ## Collapse kmers
  max_edit <- 2
  ## Find N using kmer distances 
  small_window_for_N_count <- 1000
  small_window_step_for_N_count <- 100
  small_window_min_percentage_of_distances <- small_window_for_N_count / 100

  ### Extract kmers =======================================================================================
  start_fasta_relative <- start
  end_fasta_relative <- end
  start <- 1
  end <- end_fasta_relative - start_fasta_relative + 1

  # kmers_list <- unlist(lapply(X = (start : (end - kmer)), FUN = extract_kmers, kmer, sequence))
  kmers_list <- unlist(lapply(X = (start : (end - kmer)), function(X) return(paste(sequence[X : (X + kmer - 1)], collapse = ""))))

  ### Find breaks  ========================================================================================
  window_size <- max_repeat
  if (window_size < 200) window_size <- 200
  if (window_size > end) window_size <- end

  window_starts <- NULL
  window_ends <- NULL
  window_ends_compare <- NULL
  windows_comparison_score <- NULL
  array_breaks_coordinates <- NULL

  if ((end - window_size) > start) { # if it's not, windows_comparison_score is not set and arrays are not split, so else is not needed
    window_starts <- genomic_bins_starts(start = start, end = (end - window_size), bin_size = window_step)
    if (length(window_starts) < 2) {
      window_ends <- (end - window_size)
    } else {
      window_ends <- window_starts + window_size - 1
    }

    # NEW
    window_starts <- window_starts[window_ends <= (end - window_size / 2)]
    window_ends <- window_ends[window_ends <= (end - window_size / 2)]
    window_starts_compare <- window_ends + 1
    window_ends_compare <- window_starts_compare + window_size - 1
    window_ends_compare[window_ends_compare <= end] <- end

    windows_comparison_score <- rep(0, length(window_starts))

    for (i in seq_along(window_starts_compare)) {
      kmers_window_A <- kmers_list[window_starts[i] : window_ends[i]]
      kmers_window_B <- kmers_list[window_starts_compare[i] : window_ends_compare[i]]
      windows_comparison_score[i] <- sum(kmers_window_B %in% kmers_window_A) / length(kmers_window_B)
    }
  }
  remove(window_starts, window_ends)
  gc()
  if (length(windows_comparison_score) == 0) {
    # the region is too small to compare anything, analyse it as is
    arrays <- data.frame(start = start, end = end, seqID = seqID, numID = numID, score = -4, top_N = 0, top_5_N = "", representative = "")
  } else if (sum(windows_comparison_score > min_windows_comparison_score_to_detach_array) == 0) {
    # there is no windows above the threshold to split arrays, likely a mislabeled region or a complex one or contains
    # repeats longer than the settings allow for
    arrays <- data.frame(start = start, end = end, seqID = seqID, numID = numID, score = 0, top_N = 0, top_5_N = "", representative = "")
  } else {
    #attempt splitting and save each array, allow overlaps or not
    window_event <- "new_top"
    window_event_position <- 1
    i <- which(windows_comparison_score > min_windows_comparison_score_to_detach_array)[1] + 1
    while (i <= length(window_starts_compare)) {
      while (i <= length(window_starts_compare) && window_event[length(window_event)] == "new_top") {
        if (windows_comparison_score[i] <= min_windows_comparison_score_to_split_array) {
          window_event <- c(window_event, "new_bottom")
          window_event_position <- c(window_event_position, i)
        }
        i <- i + 1
      }
      while (i <= length(window_starts_compare) && window_event[length(window_event)] == "new_bottom") {
        if (windows_comparison_score[i] >= min_windows_comparison_score_to_detach_array) {
          window_event <- c(window_event, "new_top")
          window_event_position <- c(window_event_position, i)
        }
        i <- i + 1
      }
    }

    if (length(window_event) <= 2) {
      # no second array found, report as is
      arrays <- data.frame(start = start, end = end, seqID = seqID, numID = numID, score = 0, top_N = 0, top_5_N = "", representative = "")
    } else {
      # there are multiple arrays, find them all
      array_breaks <- NULL
      for (i in 1 : (length(window_event) - 1)) {
        if (window_event[i] == "new_bottom" && window_event[i + 1] == "new_top") {
          array_breaks <- c(array_breaks, (window_event_position[i] + which.min(windows_comparison_score[window_event_position[i] : window_event_position[i + 1]]) - 1))
        }
      }
      array_breaks_coordinates <- window_ends[array_breaks] + 1
      arrays <- data.frame(start = 1,
                           end = (array_breaks_coordinates[1] - 1 + array_overlaps),
                           seqID = seqID,
                           numID = numID,
                           score = mean(windows_comparison_score[1 : array_breaks[1]]), # This score is not perfect (especially with short arrays), don't use too much
                           top_N = 0,
                           top_5_N = "",
                           representative = "")
      if (length(array_breaks) > 1) {
        for (i in 1 : (length(array_breaks) - 1)) {
          arrays <- rbind(arrays, list(start = array_breaks_coordinates[i] - array_overlaps,
                                       end = (array_breaks_coordinates[i + 1] - 1 + array_overlaps),
                                       seqID = seqID,
                                       numID = numID,
                                       score = 0,
                                       top_N = 0,
                                       top_5_N = "",
                                       representative = ""))
        }
      }
      arrays <- rbind(arrays, list(start = array_breaks_coordinates[length(array_breaks)],
                                   end = end,
                                   seqID = seqID,
                                   numID = numID,
                                   score = 0,
                                   top_N = 0,
                                   top_5_N = "",
                                   representative = ""))
    }
  }
  remove(window_starts_compare, window_ends_compare, array_breaks_coordinates)
  gc()
  if (!inherits(arrays, "data.frame")) { # sanity check, this should not happen
    arrays <- data.frame(start = start, end = end, seqID = seqID, numID = numID, score = -3, top_N = 0, top_5_N = "", representative = "")
  } else if (nrow(arrays) == 0) {
    arrays <- data.frame(start = start, end = end, seqID = seqID, numID = numID, score = -2, top_N = 0, top_5_N = "", representative = "")
  }
  ### For each array  =====================================================================================
  for (i in seq_len(nrow(arrays))) {
    arrays$score[i] <- 100 - seq_win_score_int(start = 1, end = (arrays$end[i] - arrays$start[i] + 1), kmer = kmer, fasta_extraction = sequence[arrays$start[i] : arrays$end[i]], fraction_p = 0.5)
    time_report_df <- as.numeric(arrays$end[i] - arrays$start[i])

    ## extract kmers ==========================================================================
    kmers_list_local <- kmers_list[arrays$start[i] : (arrays$end[i] - kmer)]
    counts_kmers <- table(kmers_list_local)
    counts_kmers <- counts_kmers[order(counts_kmers, decreasing = TRUE)]
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

    ## Filter kmers ===========================================================================
    time_report_df <- c(time_report_df, as.numeric(Sys.time())) # start
    kmer_names <- kmer_names[counts_kmers >= global_min_kmers_count]
    counts_kmers <- counts_kmers[counts_kmers >= global_min_kmers_count]
    if (length(counts_kmers) == 0) {
      next
    }

    min_kmers_count <- ceiling(length(counts_kmers) %/% 1000) + 1  # second filtering if there's still too many unique kmers
    if (min_kmers_count > global_min_kmers_count) {
      kmer_names_t <- kmer_names[counts_kmers >= min_kmers_count]
      counts_kmers_t <- counts_kmers[counts_kmers >= min_kmers_count]
      if (length(counts_kmers_t) == 0) {
        # Too stringent filtering, using top 10% kmers
        kmer_names <- kmer_names[1 : seq_len(ceiling(length(kmer_names) / 10))]
        counts_kmers <- counts_kmers[1 : seq_len(ceiling(length(counts_kmers) / 10))]
      } else if (length(counts_kmers_t) > 5000) {
        # Still too many unique kmers, hard cap at 5000
        kmer_names <- kmer_names[1 : 5000]
        counts_kmers <- counts_kmers[1 : 5000]
      } else {
        kmer_names <- kmer_names_t
        counts_kmers <- counts_kmers_t
      }
    }

    ## collapse kmers =========================================================================
    time_report_df <- c(time_report_df, as.numeric(Sys.time())) # fil kmers
    collapsed_kmers <- collapse_kmers(counts_kmers, kmer_names, max_edit, verbose = FALSE)
    for (j in seq_along(collapsed_kmers)) {
      collapsed_kmers[[j]]$locations <- which(kmers_list_local %in% collapsed_kmers[[j]]$kmers) - 1 + arrays$start[i]
    }
    remove(kmers_list_local, counts_kmers, kmer_names)
    gc()
    ## calculate distances ====================================================================
    time_report_df <- c(time_report_df, as.numeric(Sys.time())) # col km
    distances <- NULL
    kmer_starts <- NULL
    for (j in seq_along(collapsed_kmers)) {
      collapsed_kmers[[j]]$distances <- (collapsed_kmers[[j]]$locations[2:length(collapsed_kmers[[j]]$locations)] -
                                           collapsed_kmers[[j]]$locations[1:(length(collapsed_kmers[[j]]$locations) - 1)])
      collapsed_kmers[[j]]$locations <- collapsed_kmers[[j]]$locations[-length(collapsed_kmers[[j]]$locations)]
      collapsed_kmers[[j]]$locations <- collapsed_kmers[[j]]$locations[collapsed_kmers[[j]]$distances <= max_repeat]
      collapsed_kmers[[j]]$distances <- collapsed_kmers[[j]]$distances[collapsed_kmers[[j]]$distances <= max_repeat]
      collapsed_kmers[[j]]$locations <- collapsed_kmers[[j]]$locations[collapsed_kmers[[j]]$distances >= min_repeat]
      collapsed_kmers[[j]]$distances <- collapsed_kmers[[j]]$distances[collapsed_kmers[[j]]$distances >= min_repeat]
      kmer_starts <- c(kmer_starts, collapsed_kmers[[j]]$locations)
      distances <- c(distances, collapsed_kmers[[j]]$distances)
    }
    kmer_starts <- kmer_starts[distances <= max_repeat]
    distances <- distances[distances <= max_repeat]
    kmer_starts <- kmer_starts[distances >= min_repeat]
    distances <- distances[distances >= min_repeat]
    # add reverse distances
    kmer_starts_2 <- kmer_starts + distances

    ## Find N using kmer distances ============================================================
    time_report_df <- c(time_report_df, as.numeric(Sys.time())) # cal dist
    window_starts <- genomic_bins_starts(start = arrays$start[i], end = arrays$end[i], bin_size = small_window_step_for_N_count)

    if (length(window_starts) < 2) {
      window_ends <- arrays$end[i]
    } else {
      window_ends <- c((window_starts[2:length(window_starts)] - 1), arrays$end[i])
    }
    if (length(window_ends) != length(window_starts)) window_ends <- arrays$end[i]
    if ((window_ends[length(window_ends)] - window_starts[length(window_starts)]) < (small_window_step_for_N_count / 2)) {
      window_ends[length(window_ends) - 1] <- window_ends[length(window_ends)]
      window_ends <- window_ends[-length(window_ends)]
      window_starts <- window_starts[-length(window_starts)]
    }

    window_ends <- window_ends - small_window_step_for_N_count + small_window_for_N_count
    window_ends[window_ends > arrays$end[i]] <- arrays$end[i]

    time_report_df <- c(time_report_df, as.numeric(Sys.time())) #N A
    #TODO: optimise: This is by far the longest part, over 10% of the single array calculation and is increasing exponentially

    moving_top_distance <- vector(length = length(window_starts), mode = "numeric")
    for (j in seq_along(window_starts)) {
      which_distances <- (kmer_starts >= window_starts[j] & kmer_starts <= window_ends[j]) | (kmer_starts_2 >= window_starts[j] & kmer_starts_2 <= window_ends[j])
      if (sum(which_distances) > small_window_min_percentage_of_distances) {
        moving_top_distance[j] <- as.numeric(names(which.max(table(distances[which_distances]))))
      }
      kmer_starts <- kmer_starts[!which_distances]
      kmer_starts_2 <- kmer_starts_2[!which_distances]
      distances <- distances[!which_distances]
    }
    time_report_df <- c(time_report_df, as.numeric(Sys.time())) #N B
    remove(kmer_starts, kmer_starts_2, distances, window_starts, window_ends)
    gc()

    top_N_array <- sort(table(moving_top_distance[moving_top_distance != 0]), decreasing = TRUE)
    remove(moving_top_distance)
    gc()
    # merge N values that are only 1 bp apart,
    # TODO: consider removing low count N values in case a region contains them all,
    # or limit to top X Ns to be merged only
    # or limit to Ns with at least a min count to be allowed to be merged
    top_N_array <- top_N_array[order(as.numeric(names(top_N_array)))]
    j <- 1
    last_name <- as.numeric(names(top_N_array)[j])
    while (j < length(top_N_array)) {
      if (as.numeric(names(top_N_array)[j + 1]) == (last_name + 1)) {
        last_name <- as.numeric(names(top_N_array)[j + 1])
        top_N_array[j] <- top_N_array[j] + top_N_array[j + 1]
        names(top_N_array)[j] <- paste0(names(top_N_array)[j], ",", as.character(last_name))
        top_N_array <- top_N_array[-(j + 1)]
      } else {
        last_name <- as.numeric(names(top_N_array)[j + 1])
        j <- j + 1
      }
    }
    top_N_array <- sort(top_N_array, decreasing = TRUE)
    time_report_df <- c(time_report_df, as.numeric(Sys.time())) #N C

    top_N_distances <- NULL
    count_Ns <- length(top_N_array)

    if (count_Ns != 0) {
      if (count_Ns >= 5) count_Ns <- 5
      for (j in seq_len(count_Ns)) {
        arrays$top_5_N[i] <- paste0(arrays$top_5_N[i], "N_", names(top_N_array[j]), "_Count_", top_N_array[j], ".")
      }
      top_N_distances <- strsplit(arrays$top_5_N[i], split = "[.]")[[1]][1]
      top_N_distances <- strsplit(top_N_distances, split = "[_]")[[1]][2]
      top_N_distances <- strsplit(top_N_distances, split = "[,]")[[1]]
      top_N_distances <- as.numeric(top_N_distances)
      arrays$top_N[i] <- floor(mean(top_N_distances))
    }
    remove(top_N_array, count_Ns)
    gc()

    ## Identify kmers likely forming the repeat ===============================================
    # Use the best kmer and extract up to max_repeats_to_align, align and get consensus
    time_report_df <- c(time_report_df, as.numeric(Sys.time())) #N D
    max_repeats_to_align <- 10
    collapsed_kmers_topN_counts <- NULL
    collapsed_kmers_topN_ratio <- NULL
    for (j in seq_along(collapsed_kmers)) {
      distances_of_collapsed_kmer <- table(collapsed_kmers[[j]]$distances)
      collapsed_kmers[[j]]$count_kmers_top_N_distance <- sum(distances_of_collapsed_kmer[as.numeric(names(distances_of_collapsed_kmer)) %in% top_N_distances])
      collapsed_kmers_topN_counts <- c(collapsed_kmers_topN_counts, collapsed_kmers[[j]]$count_kmers_top_N_distance)
      collapsed_kmers_topN_ratio <- c(collapsed_kmers_topN_ratio, (collapsed_kmers[[j]]$count_kmers_top_N_distance / collapsed_kmers[[j]]$count))
    }

    top_kmer <- NULL
    top_kmer <- collapsed_kmers[[which.max(collapsed_kmers_topN_counts)]]
    top_kmer$locations <- top_kmer$locations[top_kmer$distances %in% top_N_distances]
    top_kmer$distances <- top_kmer$distances[top_kmer$distances %in% top_N_distances]

    remove(collapsed_kmers, collapsed_kmers_topN_counts, collapsed_kmers_topN_ratio, top_N_distances)
    gc()

    time_report_df <- c(time_report_df, as.numeric(Sys.time())) # identify kmers A

    # TODO: OPTIMISE? this section below is responsible for a big chink of the runtime of this function
    if(length(top_kmer$distances) > 0) {
      arrays$top_N[i] <- floor(median(top_kmer$distances))
      top_kmer_list <- NULL
      for (j in seq_along(top_kmer$locations)) {
        top_kmer_list <- c(top_kmer_list, paste(sequence[top_kmer$locations[j] : (top_kmer$locations[j] + (top_kmer$distances[j] - 1))], collapse = ""))
      }

      if (length(top_kmer_list) == 0) {
        arrays$representative[i] <- ""
      } else {
        if (length(top_kmer_list) > max_repeats_to_align) top_kmer_list <- top_kmer_list[runif(max_repeats_to_align, 1, length(top_kmer_list))]
        if (length(top_kmer_list) == 1) {
          consensus <- top_kmer_list
        } else {
          alignment <- NULL
          consensus <- NULL
          alignment <- write_align_read(mafft_exe = mafft,
                                        temp_dir = temp_dir,
                                        sequences = top_kmer_list,
                                        name = paste(seqID, numID, arrID, i, arrays$start[i], runif(1, 0, 1), sep = "_"))

          # TODO: maybe check internal duplication of the representative, to split if needed. Symmetrically (so AA into A) or assumetrically (ABB into A B and B)
          consensus <- consensus_N(alignment, arrays$top_N[i])
          remove(alignment)
          gc()
        }
        arrays$representative[i] <- consensus
        remove(top_kmer_list, consensus)
        gc()
      }
    } else {
      arrays$top_N[i] <- 0
      arrays$representative[i] <- ""
    }
    remove(top_kmer, max_repeats_to_align)
    gc()
    time_report_df <- c(time_report_df, as.numeric(Sys.time()))  # identify kmers B
    # write.csv(time_report_df, paste("Array_06_times", i, seqID, numID, arrID, ".csv", sep = "_"))
  }

  ### Prepare the output ==================================================================================
  arrays$start <- arrays$start + start_fasta_relative - 1
  arrays$end <- arrays$end + start_fasta_relative - 1
  # sink()
  # if(!sink_output) file.remove(file.path(temp_dir, paste0(seqID, "_", arrID, "_logfile.txt")))
  remove(start, end, sequence, seqID, numID, max_repeat, min_repeat, mafft, temp_dir, src_dir, kmers_list, start_fasta_relative,
  end_fasta_relative, window_step, min_windows_comparison_score_to_detach_array, min_windows_comparison_score_to_split_array, array_overlaps,
  global_min_kmers_count, max_edit, small_window_for_N_count, small_window_step_for_N_count, small_window_min_percentage_of_distances)
  gc()
  return(arrays)
}