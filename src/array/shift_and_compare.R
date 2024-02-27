shift_and_compare = function(sequence, templates = "") {
  if(sequence == "") {
    return(paste0("_", sequence))
  }
  score_threshold = 0.8

  ### Shift =================================================
  #TODO: fix this for short sequences
  sequence = shift_sequence(sequence)

  if(templates != "") {
  ### Compare ===============================================
    scores = data.frame(score = vector(mode = "numeric", length = length(templates)), 
                        sequence_shifted = vector(mode = "character", length = length(templates)))
    for (i in seq_along(templates)) {
      scores[i,] = list(compare_circular(sequence, templates[[i]])) 
    }
    if(max(scores$score) >= score_threshold) {
      sequence = scores$sequence_shifted[which.max(scores$score)]
      sequence = paste0(names(templates)[which.max(scores$score)], "_", sequence)
    }
  } else {
    sequence = paste0("_", sequence)
  }

  return(sequence)
}
