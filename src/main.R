main <- function(cmd_arguments) {
  cl <- makeCluster(cmd_arguments$cores_no)
  registerDoParallel(cl)

  # Load fasta
  fasta_content <- read_fasta_and_list(cmd_arguments$fasta_file)
  # Calculate repeat scores for windows in each sequence
  repeat_scores <- list()
  for (i in seq_along(fasta_content)) {
    append(repeat_scores, sequence_window_score(fasta_content[i], cmd_arguments$cores_no)) # nolint par_f
  }
  # Make a regions data frame
  repetitive_regions <- data.frame(start = NULL, end = NULL, seqID = NULL, score = NULL) # nolint
  for (i in seq_along(repeat_scores)) {
    regions_of_sequence <- merge_windows(repeat_scores[[i]])
    if (nrow(regions_of_sequence) != 0) {
      repetitive_regions <- rbind(repetitive_regions, regions_of_sequence)
    }
  }
  # Split regions into arrays
  arrays = foreach (i = 1 : nrow(repetitive_regions)) %dopar% { # nolint
    return(split_and_check_arrays(repetitive_regions$start[i], repetitive_regions$end[i], fasta_content[repetitive_regions$seqID[i]]))
  }
  stopCluster(cl)
}
