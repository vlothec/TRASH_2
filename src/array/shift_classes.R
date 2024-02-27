shift_classes <- function(arrays) {
  # start	end	seqID	numID	score	top_N	top_5_N	representative	class
  arrays$width <- arrays$end - arrays$start
  arrays$importance <- arrays$width * arrays$score
  if(sum(which(arrays$score == 0)) > 0) arrays$importance[arrays$score == 0] <- 0.5 * arrays$width[arrays$score == 0]
  class_sequence <- arrays$representative[which.max(arrays$importance)]
  sequences <- arrays$representative
  scores = data.frame(score = vector(mode = "numeric", length = length(sequences)), 
                      sequence_shifted = vector(mode = "character", length = length(sequences)))
  for (i in seq_along(sequences)) {
    scores[i,] = compare_circular(sequences[i], class_sequence)
  }
  return(scores$sequence_shifted)
}