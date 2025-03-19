classify_repeats <- function(repeat_df) {
  # names repeat_df:
  # start	end	seqID	numID	score	top_N	top_5_N	representative	class
  # repeat_df = read.csv("CP116284.1.fasta_arrays.csv")

  ## name rows with no sequence as "none_identified"
  repeat_df$class[repeat_df$representative == ""] <- "none_identified"
  ## make a temp df for easier handling
  which_to_classify <- which(repeat_df$class == "")
  repeat_df_temp <- repeat_df[which_to_classify, ]
  repeat_df_temp$width <- repeat_df_temp$end - repeat_df_temp$start
  repeat_df_temp$importance <- repeat_df_temp$width * repeat_df_temp$score
  repeat_df_temp$importance[repeat_df_temp$score == 0] <- 1 * repeat_df_temp$width[repeat_df_temp$score == 0]
  repeat_df_temp$rep_width = nchar(repeat_df_temp$representative)

  ## Classify; TODO: do in parallel if it takes too much time
  kmer = 7
  size_dif_to_check <- 0.15 # fraction of the checked rep size to check against
  max_distance_to_classify <- 0.85 # kmer in kmer based score threshold to classify as similar

  # Make lists of kmer vectors to help with 
  repeat_df_kmers_fw = NULL
  repeat_df_kmers_rv = NULL
  for(i in seq_len(nrow(repeat_df_temp))) {
    string_length = nchar(repeat_df_temp$representative[i])
    sequence_ext_fw = paste0(rep(repeat_df_temp$representative[i], ceiling((string_length + kmer) / string_length)), collapse = "") 
    sequence_ext_rv = paste0(rep(rev_comp_string(repeat_df_temp$representative[i]), ceiling((string_length + kmer) / string_length)), collapse = "") 
    repeat_df_kmers_fw = append(repeat_df_kmers_fw, list(unlist(lapply(seq_len(string_length), function(X) substr(sequence_ext_fw, X, (X + kmer - 1))))))
    repeat_df_kmers_rv = append(repeat_df_kmers_rv, list(unlist(lapply(seq_len(string_length), function(X) substr(sequence_ext_rv, X, (X + kmer - 1))))))
  }
  names_iterator <- 1
  while(sum(repeat_df_temp$class == "") > 0) {
    which_top <- which.max(repeat_df_temp$importance)
    string_length = nchar(repeat_df_temp$representative[which_top])
    new_class_name <- paste0(repeat_df_temp$rep_width[which_top], "_", names_iterator)
    class_bp_range <- floor(repeat_df_temp$rep_width[which_top] * (1 - size_dif_to_check)) : ceiling(repeat_df_temp$rep_width[which_top] * (1 + size_dif_to_check))
    which_to_compare <- which((repeat_df_temp$rep_width %in% class_bp_range) & (repeat_df_temp$class == ""))
    which_to_compare = which_to_compare[!(which_to_compare == which_top)]

    if(length(which_to_compare) == 0) {
      repeat_df_temp$class[which_top] <- new_class_name
      repeat_df_temp$importance[which_top] <- 0
      names_iterator = names_iterator + 1
      next
    }
    scores_fw = unlist(lapply(which_to_compare, function(X) sum(repeat_df_kmers_fw[[which_top]] %in% repeat_df_kmers_fw[[X]])))
    scores_rv = unlist(lapply(which_to_compare, function(X) sum(repeat_df_kmers_fw[[which_top]] %in% repeat_df_kmers_rv[[X]])))
    
    distances = 1 - scores_fw / length(repeat_df_kmers_fw[[which_top]])
    distances_rv = 1 - scores_rv / length(repeat_df_kmers_fw[[which_top]])

    for(i in seq_along(distances)) { #if both have some similarity, choose one
      if(distances[i] >= distances_rv[i]) {
        distances[i] = 1
      } else {
        distances_rv[i] = 1
      }
    }
    similar = which(distances <= max_distance_to_classify)
    similar_rv = which(distances_rv <= max_distance_to_classify)
 
    if (length(similar) > 0) {
      which_new_class <- which_to_compare[similar]
      repeat_df_temp$class[which_new_class] <- new_class_name
      repeat_df_temp$importance[which_new_class] <- 0
    }
    if (length(similar_rv) > 0) {
      which_new_class <- which_to_compare[similar_rv]
      repeat_df_temp$class[which_new_class] <- new_class_name
      repeat_df_temp$importance[which_new_class] <- 0
      for(k in seq_along(similar_rv)) {
        repeat_df_temp$representative[which_new_class[k]] <- rev_comp_string(repeat_df_temp$representative[which_new_class[k]])
      }
    }
    repeat_df_temp$class[which_top] <- new_class_name
    repeat_df_temp$importance[which_top] <- 0
    names_iterator = names_iterator + 1
  }

  repeat_df$class[which_to_classify] <- repeat_df_temp$class

  return(repeat_df[c("start", "end", "seqID", "numID", "score", "top_N", "top_5_N", "representative", "class")])
}