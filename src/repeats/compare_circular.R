compare_circular = function(sequence, template) {
  # brute force it with adist to all shifts
  string_length <- nchar(sequence)
  start_substrings <- lapply(seq_len(string_length), function(X) substr(sequence, 1, X))
  end_substrings <- lapply(seq_len(string_length), function(X) substr(sequence, (X + 1), string_length))
  shifted_strings <- lapply(seq_len(string_length), function(X) paste0(end_substrings[X], start_substrings[X]))
  shifted_strings_reverse <- unlist(lapply(shifted_strings, rev_comp_string))
  shifted_strings_all = c(shifted_strings, shifted_strings_reverse)
  scores = adist(template, shifted_strings_all)
  score = min(scores)
  sequence_shifted = shifted_strings_all[which.min(shifted_strings_all)]
  return(c(score, sequence_shifted))
}