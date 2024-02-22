shift_and_compare = function(sequence, templates = "") {
  if(sequence == "") {
    return(paste0("_", sequence))
  }
  score_threshold = 0.8
  assigned_to_template = FALSE

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
        assigned_to_template = TRUE
      }
  } 
  if(!assigned_to_template) {
    ### Shift =================================================
    sequence = shift_sequence(sequence) # TODO 
    sequence = paste0("_", sequence)
  } 

  if(FALSE) {
    # Do it with ffs
    if(templates != "") {
      scores = ffs_scores(sequence, templates)
      if(min(scores) <= score_threshold) {
        sequence = align_start_positions(sequence, templates[which.min(scores)]) # TODO
        sequence = paste0(names(templates)[which.min(scores)], "_", sequence)
        assigned_to_template = TRUE
      }
    }
  }
  return(sequence)
}