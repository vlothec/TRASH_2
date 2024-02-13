main <- function(cmd_arguments) {
  cl <- makeCluster(cmd_arguments$cores_no)
  registerDoParallel(cl)
  window_size <- cmd_arguments$max_rep_size * 1.2

  # Load fasta
  fasta_content <- read_fasta_and_list(cmd_arguments$fasta_file)
  # Calculate repeat scores for windows in each sequence
  repeat_scores <- list()
  for (i in seq_along(fasta_content)) {
    append(repeat_scores, sequence_window_score(fasta_content[i], window_size)) # nolint par_f
  }
  # Make a regions data frame
  repetitive_regions <- data.frame(starts = NULL, ends = NULL, scores = NULL, seqID = NULL, numID = NULL) # nolint
  for (i in seq_along(repeat_scores)) {
    regions_of_sequence <- merge_windows(repeat_scores[[i]], window_size, length(fasta_content[i]))
    if (!is.null(regions_of_sequence)) {
      regions_of_sequence$seqID <- names(fasta_content[i])
      regions_of_sequence$numID <- i
      repetitive_regions <- rbind(repetitive_regions, regions_of_sequence)
    }
  }
  # Split regions into arrays
  arrays <- foreach (i = 1 : nrow(repetitive_regions)) %dopar% { # nolint
    return(split_and_check_arrays(start = repetitive_regions$starts[i],
                                  end = repetitive_regions$ends[i],
                                  sequence = fasta_content[repetitive_regions$numID[i]][repetitive_regions$starts[i] : repetitive_regions$ends[i]]),
                                  seqID = repetitive_regions$seqID[i],
                                  numID = repetitive_regions$numID[i],
                                  window_size = (cmd_arguments$max_rep_size * 2))
  }






  stopCluster(cl)
}
# Functions TODO:

# split_and_check_arrays: identifies whether a region is a single array or multiple arrays and whether they are simple or complex, returns a list of arrays with these information, maybe N?
#