classify_repeats <- function(repeat_df) {
  # names repeat_df:
  # start	end	seqID	numID	score	top_N	top_5_N	representative	class

  ## name rows with no sequence as "none_identified"
  repeat_df$class[repeat_df$representative == ""] <- "none_identified"
  ## make a temp df for easier handling
  which_to_classify <- which(repeat_df$class == "")
  repeat_df_temp <- repeat_df[which_to_classify, ]
  repeat_df_temp$width <- repeat_df_temp$end - repeat_df_temp$start
  repeat_df_temp$importance <- repeat_df_temp$width * repeat_df_temp$score
  repeat_df_temp$importance[repeat_df_temp$score == 0] <- 1 * repeat_df_temp$width[repeat_df_temp$score == 0]
  repeat_df_temp$rep_width = nchar(repeat_df_temp$representative)

  ## Classify
  size_dif_to_check <- 0.1 # fraction of the checked rep size to check against
  max_edit_to_classify <- 0.8 # kmer in kmer based score threshold to classify as similar
  names_iterator <- 1
  while(length(which(repeat_df_temp$class == "")) > 0) {
    which_top <- which.max(repeat_df_temp$importance)
    new_class_seq <- repeat_df_temp$representative[which_top]
    new_class_name <- paste0(repeat_df_temp$rep_width[which_top], "_", names_iterator)
    class_bp_range <- floor(repeat_df_temp$rep_width[which_top] * (1 - size_dif_to_check)) : ceiling(repeat_df_temp$rep_width[which_top] * (1 + size_dif_to_check))
    which_to_compare <- which((repeat_df_temp$rep_width %in% class_bp_range) & (repeat_df_temp$class == ""))
    if(length(which_to_compare) == 0) {
      repeat_df_temp$class[which_top] <- new_class_name
      repeat_df_temp$importance[which_top] <- 0
      names_iterator = names_iterator + 1
      next
    }
    # NEW:
    kmer = 10
    string_length <- nchar(new_class_seq)
    sequence_fw <- paste0(rep(new_class_seq, ceiling((string_length + kmer) / string_length)), collapse = "") 
    kmers_fw <- unlist(lapply(seq_len(string_length), function(X) substr(sequence_fw, X, (X + kmer - 1))))
    sequence_rev = rev_comp_string(new_class_seq)
    sequence_rev <- paste0(rep(sequence_rev, ceiling((string_length + kmer) / string_length)), collapse = "") 
    kmers_rev <- unlist(lapply(seq_len(string_length), function(X) substr(sequence_rev, X, (X + kmer - 1))))
    scores_fw = table(unlist(lapply(seq_along(kmers_fw), function(X) grep(kmers_fw[X], repeat_df_temp$representative[which_to_compare]))))
    scores_rv = table(unlist(lapply(seq_along(kmers_rev), function(X) grep(kmers_rev[X], repeat_df_temp$representative[which_to_compare]))))

    not_in_scores <- which(!(seq_along(which_to_compare) %in% as.numeric(names(scores_fw))))
    for(i in seq_along(not_in_scores)) {
      scores_fw = append(scores_fw, 0)
      names(scores_fw)[length(names(scores_fw))] = not_in_scores[i]
    }
    not_in_scores <- which(!(seq_along(which_to_compare) %in% as.numeric(names(scores_rv))))
    for(i in seq_along(not_in_scores)) {
      scores_rv = append(scores_rv, 0)
      names(scores_rv)[length(names(scores_rv))] = not_in_scores[i]
    }
    scores_fw = scores_fw[order(as.numeric(names(scores_fw)))]
    scores_rv = scores_rv[order(as.numeric(names(scores_rv)))]
    distances = 1 - scores_fw / length(kmers_fw)
    distances_rv = 1 - scores_rv / length(kmers_rev)

    # OLD:
    # distances = adist(new_class_seq, repeat_df_temp$representative[which_to_compare]) / repeat_df_temp$rep_width[which_top]
    # distances_rv = adist(rev_comp_string(new_class_seq), repeat_df_temp$representative[which_to_compare]) / repeat_df_temp$rep_width[which_top]
    # max_edit_to_classify = 0.3

    similar = which(distances <= max_edit_to_classify)
    similar_rv = which(distances_rv <= max_edit_to_classify)
    if(sum(similar %in% similar_rv) > 0) {
      same_indexes <- which(similar %in% similar_rv)
      remove_fw = NULL
      remove_rv = NULL
      for(k in same_indexes) {
        if(distances[similar[k]] <= distances_rv[similar_rv[k]]) {
          remove_rv = c(remove_rv, k)
        } else {
          remove_fw = c(remove_fw, k)
        }
      }
      if (!is.null(remove_fw)) similar = similar[-remove_fw]
      if (!is.null(remove_rv)) similar_rv = similar_rv[-remove_rv]
    }

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