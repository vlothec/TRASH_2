handle_overlaps <- function(repeat_table, overlap_threshold = 0.1) {
  if (nrow(repeat_table) == 1) {
    return(repeat_table[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class")])
  }
  repeat_table <- repeat_table[order(repeat_table$start, decreasing = FALSE), ]
  repeat_table$overlap_with_next <- 0
  repeat_table$overlap_with_next[1 : (nrow(repeat_table) - 1)] <- repeat_table$end[1 : (nrow(repeat_table) - 1)] - repeat_table$start[2 : nrow(repeat_table)] + 1
  repeat_table$overlap_with_next[repeat_table$overlap_with_next < 0] <- 0

  ### First remove big overlapping
  while (sum(repeat_table$overlap_with_next) > 0) {
    i <- which.max(repeat_table$overlap_with_next)
    overlap_fraction <- repeat_table$overlap_with_next[i] / min(repeat_table$width[i : (i + 1)])
    if (overlap_fraction >= overlap_threshold) { # for big overlaps, remove one repeat
      if ((-1) %in% repeat_table$eval[i : (i + 1)]) { # default handling, lower score is worse
        remove_row = c(i, i + 1)[which.min(c(repeat_table$score[i], repeat_table$score[i + 1]))]
        repeat_table = repeat_table[-remove_row, ]
      } else { # nhmmer eval handling, higher score is worse
        remove_row = c(i, i + 1)[which.max(c(repeat_table$eval[i], repeat_table$eval[i + 1]))]
        repeat_table = repeat_table[-remove_row, ]
      }
    } else { # for smaller overlaps, divide overlap in two
      repeat_table$start[i + 1] = repeat_table$start[i + 1] + floor(repeat_table$overlap_with_next[i] / 2)
      repeat_table$score[i + 1] = 0
      repeat_table$eval[i + 1] = -1
      repeat_table$end[i] = repeat_table$end[i] - ceiling(repeat_table$overlap_with_next[i] / 2)
      repeat_table$score[i] = 0
      repeat_table$eval[i] = -1
    }
    repeat_table$overlap_with_next <- 0
    if(nrow(repeat_table) > 1) {
      repeat_table$overlap_with_next[1 : (nrow(repeat_table) - 1)] <- repeat_table$end[1 : (nrow(repeat_table) - 1)] - repeat_table$start[2 : nrow(repeat_table)] + 1
    }
    repeat_table$overlap_with_next[repeat_table$overlap_with_next < 0] <- 0

  }
  return(repeat_table[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class")])
}
# if (nrow(repeat_table) == 0) {
#   return(data.frame(seqID = vector(mode = "character"),
#                     arrayID = vector(mode = "numeric"),
#                     start = vector(mode = "numeric"),
#                     end = vector(mode = "numeric"),
#                     strand = vector(mode = "character"),
#                     score = vector(mode = "numeric"),
#                     eval = vector(mode = "numeric"),
#                     width = vector(mode = "numeric"),
#                     class = vector(mode = "character")))
# }