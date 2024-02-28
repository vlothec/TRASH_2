handle_overlaps <- function(repeat_table, overlap_threshold = 0.1, representative_len) {
  repeat_table <- repeat_table[order(repeat_table$start, decreasing = FALSE), ]
  repeat_table$overlap_with_next <- 0
  repeat_table$overlap_with_next[1 : (nrow(repeat_table) - 1)] <- repeat_table$end[1 : (nrow(repeat_table) - 1)] - repeat_table$start[2 : nrow(repeat_table)] + 1
  repeat_table$overlap_with_next[repeat_table$overlap_with_next < 0] <- 0

  while(sum(repeat_table$overlap_with_next) > 0) {
    this_smaller_than_next <- repeat_table$width[2 : nrow(repeat_table)] < repeat_table$width[1 : (nrow(repeat_table) - 1)]
    repeat_table$this_next_min_width <- repeat_table$width
    repeat_table$this_next_min_width[this_smaller_than_next] <- repeat_table$width[this_smaller_than_next] 
    repeat_table$this_next_min_width[!this_smaller_than_next] <- repeat_table$width[which(!this_smaller_than_next) + 1]

    repeat_table$overlap_fraction = repeat_table$overlap_with_next / repeat_table$this_next_min_width

    repeats_to_remove = NULL
    if(min(repeat_table$eval)[1] != -1) { # this means it was done with nhmmer and not recalculated
      while(sum(repeat_table$overlap_fraction) > 0) {
        i = which.max(repeat_table$overlap_fraction)
        if(repeat_table$overlap_fraction[i] >= overlap_threshold) {
          # mark the one with higher score for removal and 0 its score
          repeats_to_remove = c(repeats_to_remove, (which.max(c(repeat_table$eval[i], repeat_table$eval[i + 1])) - 1 + i))
          repeat_table$overlap_fraction[i] = 0
          repeat_table$overlap_fraction[i + 1] = 0
          repeat_table$score[repeats_to_remove[length(repeats_to_remove)]] = 0
        } else {
          # truncate both and adjust the score # TODO: consider truncating the one with lower score or one upstream
          repeat_table$start[i + 1] = repeat_table$start[i + 1] + floor(repeat_table$overlap_with_next[i] / 2)
          repeat_table$score[i + 1] = 0
          repeat_table$eval[i + 1] = -1
          repeat_table$end[i] = repeat_table$end[i] - ceiling(repeat_table$overlap_with_next[i] / 2)
          repeat_table$score[i] = 0
          repeat_table$eval[i] = 999
          repeat_table$overlap_fraction[i] = 0
          repeat_table$overlap_fraction[i + 1] = 0
        }
      }
    } else { # and this with default mapping or nhmmer that was recalculated
      while(sum(repeat_table$overlap_fraction) > 0) {
        i = which.max(repeat_table$overlap_fraction)
        if (repeat_table$overlap_fraction[i] >= overlap_threshold) {
          # mark the one with higher score for removal and 0 its score
          repeats_to_remove = c(repeats_to_remove, (which.max(c(repeat_table$score[i], repeat_table$score[i + 1])) - 1 + i))
          repeat_table$overlap_fraction[i] = 0
          repeat_table$overlap_fraction[i + 1] = 0
          repeat_table$score[repeats_to_remove[length(repeats_to_remove)]] = 0
        } else {
          # truncate both and adjust the score # TODO: consider truncating the one with lower score 
          repeat_table$start[i + 1] = repeat_table$start[i + 1] + floor(repeat_table$overlap_with_next[i] / 2)
          repeat_table$score[i + 1] = 0
          repeat_table$end[i] = repeat_table$end[i] - ceiling(repeat_table$overlap_with_next[i] / 2)
          repeat_table$score[i] = 0
          repeat_table$overlap_fraction[i] = 0
          repeat_table$overlap_fraction[i + 1] = 0
        }
      }
    }
    if (length(repeats_to_remove) > 0) repeat_table = repeat_table[-repeats_to_remove, ]
    if (nrow(repeat_table) == 0) {
      return(data.frame(seqID = vector(mode = "numeric"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric"),
                        eval = vector(mode = "numeric"),
                        width = vector(mode = "numeric"),
                        class = vector(mode = "character")))
    }
    repeat_table <- repeat_table[order(repeat_table$start, decreasing = FALSE), ]
    repeat_table$overlap_with_next <- 0
    if(nrow(repeat_table) > 1) {
      repeat_table$overlap_with_next[1 : (nrow(repeat_table) - 1)] <- repeat_table$end[1 : (nrow(repeat_table) - 1)] - repeat_table$start[2 : nrow(repeat_table)] + 1
    }
    repeat_table$overlap_with_next[repeat_table$overlap_with_next < 0] <- 0
  }

 
  return(repeat_table[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class")])
}