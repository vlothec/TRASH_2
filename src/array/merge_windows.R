merge_windows <- function(list_of_scores, window_size, sequence_full_length, log_messages) {
  # TODO make the treshold dynamic
  threshold <- 90

  # if(sequence_full_length < window_size) {
  #   return(data.frame(starts = vector(mode = "numeric"),
  #                     ends = vector(mode = "numeric"),
  #                     scores = vector(mode = "numeric")))
  # }

  if (sum(list_of_scores < threshold) == 0) {
    return(data.frame(starts = vector(mode = "numeric"),
                      ends = vector(mode = "numeric"),
                      scores = vector(mode = "numeric")))
  }

  # window definition needs to be propagated from sequence_window_score() funtion
  sliding_window_distance <- ceiling(window_size / 4)
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

  if (length(starts) == 1) {
    repetitive_regions <- data.frame(starts = start[list_of_scores < threshold],
                                     ends = sequence_full_length[list_of_scores < threshold],
                                     scores = list_of_scores[list_of_scores < threshold])
    return(repetitive_regions)
  }

  if(length(starts) != length(ends)) {
    stop("merge_windows starts != ends")
  }
  if(length(starts) != length(list_of_scores)) {
    stop("merge_windows starts != list_of_scores")
  }
  if(length(ends) != length(list_of_scores)) {
    stop("merge_windows ends != list_of_scores")
  }
  
  repetitive_regions <- data.frame(starts = starts[list_of_scores < threshold],
                                   ends = ends[list_of_scores < threshold],
                                   scores = list_of_scores[list_of_scores < threshold])
  
  if (nrow(repetitive_regions) < 2) return(repetitive_regions)

  i <- 1
  while (i < nrow(repetitive_regions)) {
    if ((repetitive_regions$ends[i] + 1) >= (repetitive_regions$starts[i + 1])) {
      repetitive_regions$scores[i] <- ((repetitive_regions$scores[i] * (repetitive_regions$ends[i] - repetitive_regions$starts[i])) +
                                         (repetitive_regions$scores[i + 1] * (repetitive_regions$ends[i + 1] - repetitive_regions$starts[i + 1]))) / (repetitive_regions$ends[i + 1] - repetitive_regions$starts[i + 1] + repetitive_regions$ends[i] - repetitive_regions$starts[i])
      repetitive_regions$ends[i] <- repetitive_regions$ends[i + 1]
      repetitive_regions <- repetitive_regions[-(i + 1), ]
      i <- i - 1
    }
    i <- i + 1
  }
  return(repetitive_regions)
}