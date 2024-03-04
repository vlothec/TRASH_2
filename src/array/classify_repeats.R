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
  repeat_df_temp$importance[repeat_df_temp$score == 0] <- 0.5 * repeat_df_temp$width[repeat_df_temp$score == 0]
  repeat_df_temp$rep_width = nchar(repeat_df_temp$representative)

  ## Classify
  size_dif_to_check <- 0.1 # fraction of the checked rep size to check against
  max_edit_to_classify <- 0.3 # again, fraction of the checked rep size to classify as similar
  names_iterator <- 1
  while(length(which(repeat_df_temp$class == "")) > 0) {
    which_top <- which.max(repeat_df_temp$importance)
    new_class_seq <- repeat_df_temp$representative[which_top]
    new_class_name <- paste0(repeat_df_temp$rep_width[which_top], "_", names_iterator)
    class_bp_range <- floor(repeat_df_temp$rep_width[which_top] * (1 - size_dif_to_check)) : ceiling(repeat_df_temp$rep_width[which_top] * (1 + size_dif_to_check))
    which_to_compare <- which((repeat_df_temp$rep_width %in% class_bp_range) & (repeat_df_temp$class != ""))
    if(length(which_to_compare) == 0) {
      repeat_df_temp$class[which_top] <- new_class_name
      repeat_df_temp$importance[which_top] <- 0
      names_iterator = names_iterator + 1
      next
    }
    distances = adist(new_class_seq, repeat_df_temp$representative[which_to_compare]) / repeat_df_temp$rep_width[which_top]
    similar = which(distances <= max_edit_to_classify)
    if (length(similar) > 0) {
      which_new_class <- which_to_compare[similar]
      repeat_df_temp$class[which_new_class] <- new_class_name
      repeat_df_temp$importance[which_new_class] <- 0
    }
    repeat_df_temp$class[which_top] <- new_class_name
    repeat_df_temp$importance[which_top] <- 0
    names_iterator = names_iterator + 1
  }
  repeat_df$class[which_to_classify] <- repeat_df_temp$class

  return(repeat_df[c("start", "end", "seqID", "numID", "score", "top_N", "top_5_N", "representative", "class")])
}