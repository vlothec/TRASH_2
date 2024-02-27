compare_circular = function(sequence, template) {
  # brute force it with adist to all shifts
  # TODO: Do it recursive, compare to 4 shifts, choose best, keep going
  string_length <- nchar(sequence)
  start_substrings <- lapply(seq_len(string_length), function(X) substr(sequence, 1, X))
  end_substrings <- lapply(seq_len(string_length), function(X) substr(sequence, (X + 1), string_length))
  shifted_strings <- lapply(seq_len(string_length), function(X) paste0(end_substrings[X], start_substrings[X]))
  shifted_strings_reverse <- unlist(lapply(shifted_strings, rev_comp_string))
  shifted_strings_all = unlist(c(shifted_strings, shifted_strings_reverse))
  scores = adist(template, shifted_strings_all, costs = list(ins = 1, del = 1, sub = 1))[1,]
  sequence_shifted = shifted_strings_all[which.min(scores)]
  return(list(score = min(scores), sequence_shifted = sequence_shifted))
}

# sequence = "accaaccttcttcttgcttctcaaagctttcatggtgtagccaaagtccgtatgagtctttggctttgtgtcttctaacaaggatacaattcttacgcctataagatccggttgcggtttaagttcttatactcaatcatacacatgacatcaagtcatattcgactccaaaacacta"
# template = "taaccaaccttcttcttgcttctcaaagctttcatggtgtagccaaagtccgtatgagtctttggctttgtgtcttctaacaaggatacaattcttacgcctataagatccggttgcggtttaagttcttatactcaatcatacacatgacatcaagtcatattcgactccaaaacac"
