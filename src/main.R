main <- function(cmd_arguments) {
  print("TRASH: workspace initialised")
  cl <- makeCluster(cmd_arguments$cores_no)
  registerDoParallel(cl)
  window_size <- ceiling(cmd_arguments$max_rep_size * 1.1)

  # Load fasta
  print(paste0("Loading the fasta file: ", basename(cmd_arguments$fasta_file)))
  fasta_content <- read_fasta_and_list(cmd_arguments$fasta_file)
  # Calculate repeat scores for windows in each sequence
  print("Calculating repeat scores")
  repeat_scores <- list()
  for (i in seq_along(fasta_content)) {
    repeat_scores <- append(repeat_scores, list(sequence_window_score(fasta_content[[i]], window_size))) # nolint par_f
  }
  # Make a regions data frame

  print("Identifying regions with high repeat content")
  repetitive_regions <- data.frame(starts = NULL, ends = NULL, scores = NULL, seqID = NULL, numID = NULL) # nolint
  for (i in seq_along(repeat_scores)) {
    regions_of_sequence <- merge_windows(list_of_scores = repeat_scores[[i]], window_size = window_size, sequence_full_length = length(fasta_content[[i]]))
    if (!is.null(regions_of_sequence)) {
      regions_of_sequence$seqID <- names(fasta_content)[[i]]
      regions_of_sequence$numID <- i
      repetitive_regions <- rbind(repetitive_regions, regions_of_sequence)
    }
  }
  if(!inherits(repetitive_regions, "data.frame")) stop("No regions with repeats identified")
  if(nrow(repetitive_regions) == 0) stop("No regions with repeats identified")

  write.csv(x = repetitive_regions, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_regions.csv")), row.names = FALSE)
  # Split regions into arrays
  paths <- .libPaths()
  clusterEvalQ(cl, library(doParallel))
  #arrays <- list()
  print("Identifying individual arrays with repeats")
  # function receives only a snippet of the sequence, and start and end coordinates that are relative to the whole fasta, and those should be used in return values
  arrays <- foreach (i = 1 : nrow(repetitive_regions), .combine = rbind, .export = c("split_and_check_arrays", "extract_kmers", "collapse_kmers", "genomic_bins_starts", "consensus_N", "write_align_read")) %dopar% { # nolint
    split_and_check_arrays(start = repetitive_regions$starts[i],
                                  end = repetitive_regions$ends[i],
                                  sequence = fasta_content[[repetitive_regions$numID[i]]][repetitive_regions$starts[i] : repetitive_regions$ends[i]],
                                  seqID = repetitive_regions$seqID[i],
                                  numID = repetitive_regions$numID[i],
                                  max_repeat = cmd_arguments$max_rep_size,
                                  mafft = "../dep/mafft-7.520-win64-signed/mafft-win/mafft.bat",
                                  temp_dir = cmd_arguments$output_folder,
                                  src_dir = getwd())
  }
  print("Arrays identified, saving the array table")

  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrays.csv")), row.names = FALSE)
  export_gff(annotations.data.frame = arrays,
             output = cmd_arguments$output_folder, 
             file.name = basename(cmd_arguments$fasta_file), 
             source = "TRASH",
             type = "Satellite_array", 
             seqid = 3, 
             start = 1, 
             end = 2, 
             score = 5,
             attributes = 6, 
             attribute.names = "Name=")

  print("Shifting array representative repeats and comparing to provided templates")
  stopCluster(cl)
}
# Functions TODO:

# split_and_check_arrays: identifies whether a region is a single array or multiple arrays and whether they are simple or complex, returns a list of arrays with these information, maybe N?
#