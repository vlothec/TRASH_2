classify_repeats = function(repeat_df, seq_colID, class_colID) {
  ## name rows with no sequence as "none_identified"
  repeat_df[, class_colID][repeat_df[, seq_colID] == ""] = "none_identified"


  return(repeat_df)
}