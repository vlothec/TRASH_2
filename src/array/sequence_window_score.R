sequence_window_score <- function(fasta_sequence, window_size, kmer = 10, output_dir = ".") {
  # TODO check if changing these settings below can make the script work better,
  # although these were optimised
  fraction_p <- 0.5
  if ((window_size / 2) <= kmer) stop("sequence_window_score: window size is too small")

  sequence_full_length <- length(fasta_sequence)

  # window definition needs to be propagated into merge_windows() funtion
  sliding_window_distance <- ceiling(window_size / 4) #if 1, it's off, use 4
  if(sliding_window_distance >= sequence_full_length) {
    starts <- 1
    ends <- sequence_full_length
  } else {
    starts <- genomic_bins_starts(start = 1, end = sequence_full_length, bin_size = sliding_window_distance)
    starts <- starts[starts < sequence_full_length]
    if (length(starts) == 1) {
      ends <- sequence_full_length
    } else {
      ends <- c((starts[2 : length(starts)] - 1), sequence_full_length) + window_size # This makes overlapping windows!
    }
    ends[ends > sequence_full_length] <- sequence_full_length
  }
  # cat(seq_along(starts), "\n")
  #divide into chunks
  scores <- NULL
  wins_per_chunk <- 100
  # make chunks starting at every 100th window start position (so default 2000 bp overlapping by 1500 bp, making around 50 Kbp), so max 100 CPUs is used in parallel
  chunk_starts <- seq(1, length(starts), wins_per_chunk)
  chunk_starts <- c(chunk_starts, (length(starts) + 1))
  cat("Sequence full length: ", round(sequence_full_length / 1000000, 3), " Mbp \t", sep = "")
  cat("Chunks to complete: ", (length(chunk_starts) - 1), ". Finished: ", sep = "")
  date <- Sys.Date()
  for(i in 1 : (length(chunk_starts) - 1)) {
    sequence_substring <- fasta_sequence[starts[chunk_starts[i]] : ends[chunk_starts[i + 1] - 1]]
    # scores <- c(scores, foreach::foreach (j = (chunk_starts[i] : (chunk_starts[i + 1] - 1)),
    foreach::foreach (j = (chunk_starts[i] : (chunk_starts[i + 1] - 1)),
                      # .combine = c,
                      .export = c("seq_win_score_int")) %dopar% {
      # cat(paste("A", i, j, "",  sep = " "))
      result <- seq_win_score_int(1, window_size, kmer, sequence_substring[(starts[j] - starts[chunk_starts[i]] + 1) : (ends[j] - starts[chunk_starts[i]] + 1)], fraction_p)
      save(result, file = paste0(output_dir, "/", i, "_", j, "_", date, "_sequence_window_score_data"))
      remove(result)
      gc()
      # return(result)
    }
    cat(i, "")
    for(j in (chunk_starts[i] : (chunk_starts[i + 1] - 1))) {
      load(paste0(output_dir, "/", i, "_", j, "_", date, "_sequence_window_score_data"))
      unlink(paste0(output_dir, "/", i, "_", j, "_", date, "_sequence_window_score_data"))
      scores <- c(scores, result)
      remove(result)
    }
  }

  cat("\n")
  # cat("a ", length(starts), sequence_full_length, "\n")
  # scores <- foreach::foreach (i = seq_along(starts),
  #                             .combine = c,
  #                             .export = c("seq_win_score_int")) %dopar% {
  #   cat(paste("A", i, "",  sep = " "))
  #   result <- seq_win_score_int(1, window_size, kmer, fasta_sequence[starts[i] : ends[i]], fraction_p)
  #   cat(paste("B", i, "",  sep = " "))
  #   return(result)
  # }
  # cat("C ")

  if ((sum(is.na(scores))>0) || (length(scores) == 0)) {
    stop("sequence_window_score did not produce valid result")
  }
  return(scores)
}
