shift_sequence <- function(sequence, k = 6) {

  string_length <- nchar(sequence)
  # extend the sequence to easily extract kmers
  sequence_fw <- paste0(rep(sequence, ceiling((string_length + k) / string_length)), collapse = "") 
  # extract the kmers
  kmers <- unlist(lapply(seq_len(string_length), function(X) substr(sequence_fw, X, (X + k - 1))))
  # calculate unique kmer hash
  # TODO: prepare a matrix that for all kmer gives a score, to not re-calculate it
  kmers_scores <- unlist(lapply(kmers, kmer_hash_score))
  kmers_scores <- c(kmers_scores, kmers_scores)
  shifts_scores_fw <- unlist(lapply(X = seq_len(string_length), function(X) sum(kmers_scores[X : (X + string_length - 1)] * seq_len(string_length))))

  # Same for the reverse comp
  sequence_rev = rev_comp_string(sequence)
  sequence_rev <- paste0(rep(sequence_rev, ceiling((string_length + k) / string_length)), collapse = "") 
  kmers <- unlist(lapply(seq_len(string_length), function(X) substr(sequence_rev, X, (X + k - 1))))
  kmers_scores <- unlist(lapply(kmers, kmer_hash_score))
  kmers_scores <- c(kmers_scores, kmers_scores)
  shifts_scores_rev <- unlist(lapply(X = seq_len(string_length), function(X) sum(kmers_scores[X : (X + string_length - 1)] * seq_len(string_length))))

 # Find the min score shift an reassign
  shift <- which.min(c(shifts_scores_fw, shifts_scores_rev))
  if(shift > string_length) {
    sequence = rev_comp_string(sequence)
    shift = shift - string_length
  }
  if(shift != 1) {
    sequence = strsplit(sequence, split = "")[[1]]
    sequence = paste(c(sequence[shift : length(sequence)], sequence[1 : (shift - 1)]), collapse = "")
  }
  return(sequence)
}