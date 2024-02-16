split_and_check_arrays <- function(start, end, sequence, seqID, numID, max_repeat, mafft, temp_dir, src_dir) {
  .libPaths(c(.libPaths(), gsub("src", "R_libs", getwd())))
  #sink(file.path(temp_dir, paste0(seqID, "_", start, "_logfile.txt")))
  print("started")
  print(start)
  print(end)
  print(max_repeat)
  ### Extract kmers and calculate distances ===============================================================
  kmer <- 12
  min_kmers_count <- 2
  max_edit <- 3

  start_fasta_relative = start
  end_fasta_relative = end
  start = 1
  end = end_fasta_relative - start_fasta_relative + 1

  print(paste0("Main kmer list of a sequence length ", length(sequence)))
  kmers_list <- unlist(lapply(X = (start : (end - kmer)), FUN = extract_kmers, kmer, sequence))
  counts_kmers <- table(kmers_list)
  kmer_names <- names(counts_kmers)

  #clean up N containing kmers
  counts_kmers <- counts_kmers[!grepl("N", kmer_names)]
  kmer_names <- kmer_names[!grepl("N", kmer_names)]
  counts_kmers <- counts_kmers[!grepl("n", kmer_names)]
  kmer_names <- kmer_names[!grepl("n", kmer_names)]

  if (length(counts_kmers) == 0) {
    return(data.frame())
  }

  kmer_names <- kmer_names[counts_kmers >= min_kmers_count]
  counts_kmers <- counts_kmers[counts_kmers >= min_kmers_count]

  if (length(counts_kmers) == 0) {
    return(data.frame())
  }

  print("Collapse kmers")
  collapsed_kmers <- collapse_kmers(counts_kmers, kmer_names, max_edit, verbose = FALSE)

  print("Find locations")
  for (i in seq_along(collapsed_kmers)) {
    collapsed_kmers[[i]]$locations <- which(kmers_list %in% collapsed_kmers[[i]]$kmers)
  }

  print("Distances")
  distances <- NULL
  kmer_starts <- NULL # these coordinates are relative to start variable
  for (i in seq_along(collapsed_kmers)) {
    collapsed_kmers[[i]]$distances <- (collapsed_kmers[[i]]$locations[2 : length(collapsed_kmers[[i]]$locations)] -
                                         collapsed_kmers[[i]]$locations[1 : (length(collapsed_kmers[[i]]$locations)-1)])
    collapsed_kmers[[i]]$locations <- collapsed_kmers[[i]]$locations[-length(collapsed_kmers[[i]]$locations)]
    collapsed_kmers[[i]]$locations <- collapsed_kmers[[i]]$locations[collapsed_kmers[[i]]$distances <= max_repeat]
    collapsed_kmers[[i]]$distances <- collapsed_kmers[[i]]$distances[collapsed_kmers[[i]]$distances <= max_repeat]
    kmer_starts <- c(kmer_starts, collapsed_kmers[[i]]$locations)
    distances <- c(distances, collapsed_kmers[[i]]$distances)
  }
  kmer_starts <- kmer_starts[distances <= max_repeat]
  distances <- distances[distances <= max_repeat]

  ### =====================================================================================================

  ### Find breaks  ========================================================================================
  print("Breaks")
  window_size <- max_repeat
  if (window_size < 500) window_size <- 500
  window_step <- 50
  min_windows_comparison_score_to_detach_array <- 0.08
  min_windows_comparison_score_to_split_array <- 0.04
  array_overlaps <- ceiling(max_repeat / 2)
  array_overlaps <- 0

  window_starts <- NULL
  window_ends <- NULL
  window_ends_compare <- NULL
  windows_comparison_score <- NULL
  array_breaks_coordinates <- NULL

  if ((end - max_repeat) > start) { # if it's not, windows_comparison_score is not set and arrays are not split, so else is not needed
    window_starts <- genomic_bins_starts(start = start, end = (end - window_size), bin_size = window_step)
    if (length(window_starts) < 2) {
      window_ends <- (end - window_size)
    } else {
      window_ends <- c((window_starts[2:length(window_starts)] - 1), (end - window_size))
    }
    if (length(window_ends) != length(window_starts)) window_ends <- (end - window_size)
    window_ends <- window_ends + window_size - window_step
    window_ends[window_ends > end] <- end

    window_starts_compare <- window_ends + 1

    if (length(window_starts_compare) < 2) {
      window_ends_compare <- end
    } else {
      window_ends_compare <- c((window_starts_compare[2:length(window_starts_compare)] - 1), end)
    }
    if (length(window_ends_compare) != length(window_starts_compare)) window_ends_compare <- end
    window_ends_compare <- window_ends_compare + window_size - window_step
    for (i in which(window_ends_compare > end)) {
      window_starts[i] <- window_starts[i] + (window_ends_compare[i] - end)

      window_ends_compare[i] <- end
    }

    windows_comparison_score <- rep(0, length(window_starts))

    for (i in seq_along(window_starts_compare)) {
      kmers_window_A <- kmers_list[window_starts[i] : window_ends[i]]
      kmers_window_B <- kmers_list[window_starts_compare[i] : window_ends_compare[i]]
      windows_comparison_score[i] <- sum(kmers_window_A %in% kmers_window_B) / window_size
    }
  }
  print("windows_comparison_score")
  print(windows_comparison_score)

  if (length(windows_comparison_score) == 0) {
    # the region is too small to compare anything, analyse it as is
    arrays <- data.frame(start = start, end = end, seqID = seqID, numID = numID, score = -4, top_N = 0, top_5_N = "", representative = "")
  } else if (sum(windows_comparison_score > min_windows_comparison_score_to_detach_array) == 0) {
    # there is no windows above the threshold to split arrays, likely a mislabeled region or a complex one or contains
    # repeats longer than the settings allow for
    arrays <- data.frame(start = start, end = end, seqID = seqID, numID = numID, score = mean(windows_comparison_score), top_N = 0, top_5_N = "", representative = "")
  } else {
    #attempt splitting and save each array, ALLOW overlaps
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
      print("window_events")
      print(window_event)
      print(window_event_position)
    
    if (length(window_event) <= 2) {
      # no second array found, report as is
      arrays <- data.frame(start = start, end = end, seqID = seqID, numID = numID, score = mean(windows_comparison_score), top_N = 0, top_5_N = "", representative = "")
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
                           score = mean(windows_comparison_score[1 : array_breaks[1]]),# This score is not perfect (especially with short arrays), don't use too much
                           top_N = 0,
                           top_5_N = "",
                           representative = "")
      if (length(array_breaks) > 1) {
        for (i in 1 : (length(array_breaks) - 1)) {
          arrays <- rbind(arrays, list(start = array_breaks_coordinates[i] - array_overlaps,
                                       end = (array_breaks_coordinates[i + 1] - 1 + array_overlaps),
                                       seqID = seqID,
                                       numID = numID,
                                       score = mean(windows_comparison_score[array_breaks[i] : array_breaks[i + 1]]),
                                       top_N = 0,
                                       top_5_N = "",
                                       representative = ""))
        }
      }
      arrays <- rbind(arrays, list(start = array_breaks_coordinates[length(array_breaks)],
                                   end = end,
                                   seqID = seqID,
                                   numID = numID,
                                   score = mean(windows_comparison_score[array_breaks[length(array_breaks)] : length(windows_comparison_score)]),
                                   top_N = 0,
                                   top_5_N = "",
                                   representative = ""))
    }
  }
  print("Arrays done")
  if (!inherits(arrays, "data.frame")) { # sanity check, this should not happen
    arrays <- data.frame(start = start, end = end, seqID = seqID, numID = numID, score = -3, top_N = 0, top_5_N = "", representative = "")
    print("Check array identification")
  } else if (nrow(arrays) == 0) {
    arrays <- data.frame(start = start, end = end, seqID = seqID, numID = numID, score = -2, top_N = 0, top_5_N = "", representative = "")
    print("Check array identification")
  }
  ### =====================================================================================================

  ### For each array  =====================================================================================
  arrays$top_N <- 0
  arrays$top_5_N <- "."
  arrays$representative <- ""
  for (i in seq_len(nrow(arrays))) {
    print(paste0("Array ", i, " / ", nrow(arrays)))
    ## recalculate distances if more than 1 array =============================================
    if (nrow(arrays) != 1) {

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

      kmer_names <- kmer_names[counts_kmers >= min_kmers_count]
      counts_kmers <- counts_kmers[counts_kmers >= min_kmers_count]
      if (length(counts_kmers) == 0) {
        print("Empty array, should not happen")
        next
      }

      collapsed_kmers <- collapse_kmers(counts_kmers, kmer_names, max_edit)
      for (j in seq_along(collapsed_kmers)) {
        collapsed_kmers[[j]]$locations <- which(kmers_list %in% collapsed_kmers[[j]]$kmers) - 1 + arrays$start[i]
      }

      distances <- NULL
      kmer_starts <- NULL
      for (j in seq_along(collapsed_kmers)) {
        collapsed_kmers[[j]]$distances <- (collapsed_kmers[[j]]$locations[2:length(collapsed_kmers[[j]]$locations)] -
                                             collapsed_kmers[[j]]$locations[1:(length(collapsed_kmers[[j]]$locations) - 1)])
        collapsed_kmers[[j]]$locations <- collapsed_kmers[[j]]$locations[-length(collapsed_kmers[[j]]$locations)]
        collapsed_kmers[[j]]$locations <- collapsed_kmers[[j]]$locations[collapsed_kmers[[j]]$distances <= max_repeat]
        collapsed_kmers[[j]]$distances <- collapsed_kmers[[j]]$distances[collapsed_kmers[[j]]$distances <= max_repeat]
        kmer_starts <- c(kmer_starts, collapsed_kmers[[j]]$locations)
        distances <- c(distances, collapsed_kmers[[j]]$distances)
      }
      kmer_starts <- kmer_starts[distances <= max_repeat]
      distances <- distances[distances <= max_repeat]
    }
    ## ========================================================================================

    ## Find N using kmer distances ============================================================
    small_window_for_N_count <- 1000
    small_window_step_for_N_count <- 100
    small_window_min_percentage_of_distances <- small_window_for_N_count / 10
    #window_starts <- genomic_bins_starts(start = arrays$start[i], end = arrays$end[i], bin_size = small_window_step_for_N_count)
    window_starts <- genomic_bins_starts(start = arrays$start[i], end = arrays$end[i], bin_size = small_window_step_for_N_count)

    if (length(window_starts) < 2) {
      window_ends <- arrays$end[i]
    } else {
      window_ends <- c((window_starts[2:length(window_starts)] - 1), arrays$end[i])
    }
    if (length(window_ends) != length(window_starts)) window_ends <- arrays$end[i]
    window_ends <- window_ends - small_window_step_for_N_count + small_window_for_N_count
    window_ends[window_ends > arrays$end[i]] <- arrays$end[i]

    moving_top_distance <- vector(length = length(window_starts), mode = "numeric")

    for (j in seq_along(window_starts)) {
      which_distances <- (kmer_starts >= window_starts[j] & kmer_starts <= window_ends[j])
      if (sum(which_distances) > small_window_min_percentage_of_distances) {
        moving_top_distance[j] <- as.numeric(names(which.max(table(distances[which_distances]))))
      }
    }
    print("moving_top_distance: ")
    print(moving_top_distance)
    top_N_array <- sort(table(moving_top_distance[moving_top_distance != 0]), decreasing = TRUE)

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

    print("top_N_array")
    print(top_N_array)

    top_N_distances <- NULL
    count_Ns <- length(top_N_array)

    if (count_Ns != 0) {
      if (count_Ns >= 5) count_Ns <- 5
      for (j in seq_len(count_Ns)) {
        arrays$top_5_N[i] <- paste0(arrays$top_5_N[i], "N_", names(top_N_array[j]), "_Count_", top_N_array[j], ".")
      }
      top_N_distances <- strsplit(arrays$top_5_N[i], split = "[.]")[[1]][2]
      top_N_distances <- strsplit(top_N_distances, split = "[_]")[[1]][2]
      top_N_distances <- strsplit(top_N_distances, split = "[,]")[[1]]
      top_N_distances <- as.numeric(top_N_distances)
      arrays$top_N[i] <- floor(mean(top_N_distances))
    }
    ## ========================================================================================

    print("arrays")
    print(arrays)
    
    ## Identify kmers likely forming the repeat ===============================================
    # Use the best kmer and extract up to 20 repeats, align and get consensus
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

    print("top_kmer")
    print(top_kmer)

    top_kmer_list <- NULL
    for (j in seq_along(top_kmer$locations)) {
      top_kmer_list <- c(top_kmer_list, paste(sequence[top_kmer$locations[j] : (top_kmer$locations[j] + (top_kmer$distances[j] - 1))], collapse = ""))
    }
    print("top_kmer_list")
    print(top_kmer_list)
    
    if (length(top_kmer_list) == 0) {
      print("No top kmers, check if as intended")
      arrays$representative[i] = ""
    } else {
      if (length(top_kmer_list) > 20) top_kmer_list <- top_kmer_list[runif(20, 1, length(top_kmer_list))]

      alignment <- NULL
      consensus <- NULL
      alignment <- write_align_read(mafft_exe = mafft,
                                    temp_dir = temp_dir,
                                    sequences = top_kmer_list,
                                    name = paste(seqID, numID, i, sep = "_"))
      consensus <- consensus_N(alignment, arrays$top_N[i])

      arrays$representative[i] <- consensus
    }
  print("Out?")
    
    ## ========================================================================================

    ## Plot arrays ============================================================================
    # if(F) {
    #   gc_calc_window = 1000

    #   par(mfrow = c(3,1), mar = c(4,4,1,1))
    #   # plot 1: moving top distance (N)
    #   plot(x = window_starts[moving_top_distance != 0], y = moving_top_distance[moving_top_distance != 0], pch = 16, 
    #        ylim = c(0, max_repeat),
    #        xlim = c(start, end))
    #   if(length(array_breaks_coordinates) != 0) abline(v = array_breaks_coordinates)
    #   # plot 2: windows comparison score
    #   plot(x = window_starts_compare[window_starts_compare > arrays$start[i] & window_starts_compare < arrays$end[i]], 
    #        y = windows_comparison_score[window_starts_compare > arrays$start[i] & window_starts_compare < arrays$end[i]], pch = 16, 
    #        ylim = c(0, 1),
    #        xlim = c(start, end))
    #   if(length(array_breaks_coordinates) != 0) abline(v = array_breaks_coordinates)
    #   # plot 3: GC
    #   window_starts <- genomic_bins_starts(start = arrays$start[i], end = arrays$end[i], bin_size = gc_calc_window)
    #   gc_track = calculate_GC_in_windows(windows.starts = (window_starts - window_starts[1] + 1), sequence = sequence[arrays$start[i] : arrays$end[i]], bin.size = gc_calc_window)
    #   plot(x = window_starts, y = gc_track, ylim = c(0,1), pch = 16, xlim = c(start, end))
    #   if(length(array_breaks_coordinates) != 0) abline(v = array_breaks_coordinates)
    #   par(mfrow = c(1,1))
    # }
    ## ========================================================================================
  }
  ### =====================================================================================================
  arrays$start = arrays$start + start_fasta_relative - 1
  arrays$end = arrays$end + start_fasta_relative - 1
  print("Returning arrays")
  print(arrays)
  return(arrays)
}
#export.gff(annotations.data.frame = arrays,output = temp_dir, file.name = "arrays.gff", seqid = "Rbrevi_chr1_extr1_edited", start = 1, end = 2, attributes = 6, attribute.names = "Length:")
#Test
#test_fasta = read_fasta_and_list("../../testing_fastas/SUPER_15_extraction.fasta")
#test_fasta = read_fasta_and_list("../../testing_fastas/ath_Chr1_extraction.fasta")
#test_fasta = read_fasta_and_list("../../testing_fastas/Rbrevi_chr1_extr1_edited.fasta")
#seqID = "sequence"
#numID = 111
#start = 1
#sequence = test_fasta[[1]]
#end = length(sequence)
#window_size = 500
#max_repeat = 400
#mafft = "C:/Users/vlothec/Documents/GitHub/TRASH_dev/temp/mafft-7.520-win64-signed/mafft-win/mafft.bat"
#temp_dir = "../../temp"
#split_and_check_arrays(1,10000,test_fasta[[1]], "CP068268", 1)
