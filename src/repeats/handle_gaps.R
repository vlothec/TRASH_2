handle_gaps <- function(repeat_table, representative_len) {
  repeat_table <- repeat_table[order(repeat_table$start, decreasing = FALSE), ]
  repeat_table$gap_to_next <- 0
  repeat_table$gap_to_next[1 : (nrow(repeat_table) - 1)] <- repeat_table$start[2 : nrow(repeat_table)] - repeat_table$end[1 : (nrow(repeat_table) - 1)] - 1
  repeat_table$gap_to_next[repeat_table$gap_to_next < 0] <- 0
  repeat_table$gap_to_next[nrow(repeat_table)] <- 9999999

  repeats_to_remove = NULL

  for(i in seq_len(nrow(repeat_table))) {
    if(repeat_table$gap_to_next[i] > 0) {
      if(repeat_table$gap_to_next[i] < (representative_len / 2)) {
        ## Extend the upstream repeat
        if(repeat_table$strand[i] == "+") {
          repeat_table$end[i] = repeat_table$end[i] + repeat_table$gap_to_next[i]
        } else {
          repeat_table$start[i + 1] = repeat_table$start[i + 1] - repeat_table$gap_to_next[i]
        }
        repeat_table$gap_to_next[i] = 0
      } else {
        ## If the element has no neighbours, remove (is either first or the previous one had a gap that's not fixed)
        if(i == 1 || repeat_table$gap_to_next[i - 1] != 0) {
          repeats_to_remove = c(repeats_to_remove, i)
        }
      }
    }
  }
  if (length(repeats_to_remove) != 0) repeat_table = repeat_table[-repeats_to_remove,]
  if (nrow(repeat_table) == 0) {
    return(data.frame(seqID = vector(mode = "character"),
                      arrayID = vector(mode = "numeric"),
                      start = vector(mode = "numeric"),
                      end = vector(mode = "numeric"),
                      strand = vector(mode = "character"),
                      score = vector(mode = "numeric"),
                      eval = vector(mode = "numeric"),
                      width = vector(mode = "numeric"),
                      class = vector(mode = "character")))
  }

  repeat_table$gap_to_next[nrow(repeat_table)] <- 0
  for(i in seq_len(nrow(repeat_table))) {
    if(repeat_table$gap_to_next[i] > 0) {
        ## Annotate the remaining gaps as an "<Repeat_name>_interspersed_element"
        repeat_table[(nrow(repeat_table) + 1), ] = list(seqID = repeat_table$seqID[i], 
                                                        arrayID = repeat_table$arrayID[i], 
                                                        start = (repeat_table$end[i] + 1), 
                                                        end = (repeat_table$start[i + 1] - 1), 
                                                        strand = ".", 
                                                        score = 100, 
                                                        eval = 0, 
                                                        width = repeat_table$gap_to_next[i],
                                                        class = paste0(repeat_table$class[i], "_interspersed_element"),
                                                        gap_to_next = 0)
    }
  }

  return(repeat_table[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class")])
}
