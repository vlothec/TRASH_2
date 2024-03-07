shift_and_compare = function(sequence, templates = 0) {
  if(sequence == "") {
    return(paste0("_", sequence))
  }
  score_threshold = 0.2
  max_size_dif = 0.15

  ### Shift =================================================
  #TODO: fix this for short sequences, for now it's set so it's not working for short
  sequence = shift_sequence(sequence)
  if(!inherits(sequence, "character")) {
    print(paste0("shift_and_compare: sequence is not character"))
    return("_")
  }

  if(!inherits(templates, "numeric")) {
  ### Compare ===============================================
    scores = data.frame(score = vector(mode = "numeric", length = length(templates)), 
                        sequence_shifted = vector(mode = "character", length = length(templates)))
    for (i in seq_along(templates)) {
      scores[i,] = compare_circular(sequence, paste(templates[[i]], collapse = ""), max_size_dif)
    }
    if(min(scores$score) <= score_threshold) {
      sequence = scores$sequence_shifted[which.min(scores$score)]
      sequence = paste0(names(templates)[which.min(scores$score)], "_", sequence)
    } else {
      sequence <- paste0("_", sequence)
    }
  } else {
    sequence = paste0("_", sequence)
  }
  if(!inherits(sequence, "character")) print(paste0("shift_and_compare: sequence is not character"))
  return(sequence)
}
