merge_windows <- function(list_of_scores, window_size, sequence_full_length) {
  # TODO make the treshold dynamic
  threshold <- 99

  if (length(list_of_scores < threshold) == 0) {
    return(NULL)
  }

  # This window definition needs to be propagated from sequence_window_score() funtion
  starts <- genomic_bins_starts(start = 1, end = sequence_full_length, bin_size = window_size)
  if (length(starts) < 2) {
    ends <- sequence_full_length
  } else {
    ends <- c((starts[2:length(starts)] - 1), sequence_full_length)
  }
  if (length(ends) != length(starts)) ends <- sequence_full_length
  if ((ends[length(ends)] - starts[length(starts)]) < (window_size / 2)) {
    ends[length(ends) - 1] <- ends[length(ends)]
    ends <- ends[-length(ends)]
    starts <- starts[-length(starts)]
  }

  repetitive_regions <- data.frame(starts = starts[list_of_scores < threshold],
                                   ends = ends[list_of_scores < threshold],
                                   scores = list_of_scores[list_of_scores < threshold])
  
  if (nrow(repetitive_regions) == 0) return(NULL)

  i <- 1
  while (i < nrow(repetitive_regions)) {
    if ((repetitive_regions$ends[i] + 1) == (repetitive_regions$starts[i + 1])) {
      repetitive_regions$scores[i] <- ((repetitive_regions$scores[i] * (repetitive_regions$ends[i] - repetitive_regions$starts[i])) +
                                         (repetitive_regions$scores[i + 1] * (repetitive_regions$ends[i + 1] - repetitive_regions$starts[i + 1]))) / (repetitive_regions$ends[i + 1] - repetitive_regions$starts[i])
      repetitive_regions$ends[i] <- repetitive_regions$ends[i + 1]
      repetitive_regions <- repetitive_regions[-(i + 1), ]
      i <- i - 1
    }
    i <- i + 1
  }
  return(repetitive_regions)
}