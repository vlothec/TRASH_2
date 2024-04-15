find_edge_best_start_end <- function(sequence_potential, representative, edge, max_start = 1) { # two sequences are character vectors

  if (edge == "right") {
    if (length(sequence_potential) < 2) return(1) #default case for right edge, returning 1 will discard new repeat in the higher function
    if (length(representative) < 2) return(1)
    new_end <- length(sequence_potential)

    scores <- rep(100, length(sequence_potential))

    # check scores in groups of 10, recursive if sequence is long (divides once length 10+, twice 100+, thrice 1000+ etc)
    scores_to_check <- round(seq(1, length(sequence_potential), length.out = 10))
    while (length(unique(scores_to_check)) == 10) {
      scores[scores_to_check] <- adist() # TODO: adist from partial potential seq to whole representative, div by partial seq length
      if (mean(scores_to_check[2 : 10] - scores_to_check[1 : 9]) > 10) {
        if (min(scores[scores_to_check]) == 1) {
          scores_to_check <- round(seq(1, scores_to_check[2], length.out = 10))
        } else if (min(scores[scores_to_check]) == length(sequence_potential)) {
          scores_to_check <- round(seq(scores_to_check[length(sequence_potential)], length(sequence_potential), length.out = 10))
        } else {
          scores_to_check <- round(seq(scores_to_check[which.min(scores[scores_to_check]) - 1], scores_to_check[which.min(scores[scores_to_check]) + 1], length.out = 10))
        }
      }
    }
    if (min(scores) != 100) {
      if (which.min(scores) == 1) {
        potential_scores <- 1 : 20
      } else if (which.min(scores) == length(sequence_potential)) {
        potential_scores <- (length(sequence_potential) - 20) : length(sequence_potential)
      } else {
        potential_scores <- (which.min(scores) - 10) : (which.min(scores) + 10)
      }
    } else {
      potential_scores <- seq_along(sequence_potential)
    }
    potential_scores <- potential_scores[potential_scores > 0]
    potential_scores <- potential_scores[potential_scores <= length(sequence_potential)]

    scores[potential_scores] <- adist() # TODO

    new_end <- which.min(scores)

    return(new_end)
  } else if (edge == "left") {
    if (length(sequence_potential) < 2) return(max_start) # default case when the repeat would be 1 bp long, returning it will discard as in "right" case above
    if (length(representative) < 2) return(max_start)
    new_start <- 1



    return(new_start)
  } else {
    stop("find_edge_best_start_end: wrond edge provided")
  }
}
