main <- function(cmd_arguments) {
  print("TRASH: workspace initialised")
  cl <- makeCluster(cmd_arguments$cores_no)
  registerDoParallel(cl)
  window_size <- ceiling(cmd_arguments$max_rep_size * 1.1)

  # Load fasta
  print(paste0("Loading the fasta file: ", basename(cmd_arguments$fasta_file)))
  fasta_content <- read_fasta_and_list(cmd_arguments$fasta_file)
  gc()
  # Calculate repeat scores for windows in each sequence
  print("Calculating repeat scores")
  repeat_scores <- list()
  for (i in seq_along(fasta_content)) {
    repeat_scores <- append(repeat_scores, list(sequence_window_score(fasta_content[[i]], window_size))) # nolint par_f
  }
  gc()
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
  gc()
  if(!inherits(repetitive_regions, "data.frame")) stop("No regions with repeats identified")
  if(nrow(repetitive_regions) == 0) stop("No regions with repeats identified")

  # write.csv(x = repetitive_regions, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_regions.csv")), row.names = FALSE)
  
  # Split regions into arrays
  clusterEvalQ(cl, library(doParallel))
  print("Identifying individual arrays with repeats")
  # function receives only a snippet of the sequence, and start and end coordinates that are relative to the whole fasta, and those should be used in return values
  arrays <- foreach (i = seq_len(nrow(repetitive_regions)), .combine = rbind, .export = c("split_and_check_arrays", "extract_kmers", "collapse_kmers", "genomic_bins_starts", "consensus_N", "write_align_read")) %dopar% {
    split_and_check_arrays(start = repetitive_regions$starts[i],
                                  end = repetitive_regions$ends[i],
                                  sequence = fasta_content[[repetitive_regions$numID[i]]][repetitive_regions$starts[i] : repetitive_regions$ends[i]],
                                  seqID = repetitive_regions$seqID[i],
                                  numID = repetitive_regions$numID[i],
                                  max_repeat = cmd_arguments$max_rep_size,
                                  min_repeat = cmd_arguments$min_rep_size,
                                  mafft = "../dep/mafft-7.520-win64-signed/mafft-win/mafft.bat",
                                  temp_dir = cmd_arguments$output_folder,
                                  src_dir = getwd(),
                                  sink_output = FALSE)
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

  arrays$representative <- foreach (i = seq_len(nrow(arrays)), .combine = c, .export = "shift_and_compare") %dopar% {
    shift_and_compare(arrays$representative[i], cmd_arguments$templates) #TODO
  }

  print("Mapping the array representative to the array using nhmmer")

  repeats <- foreach (i = seq_len(nrow(arrays)), .combine = rbind, .export = c("read_and_format_nhmmer", "handle_overlaps", "export_gff")) %dopar% {
    sink(file.path(cmd_arguments$output_folder, paste0("array_", i, "_logfile.txt")))
    print("deb0")
    .libPaths(c(.libPaths(), gsub("src", "R_libs", getwd()))) #TODO needed?
    if(arrays$representative[i] == "") {
      sink()
      gc()
      return(data.frame(seqID = vector(mode = "numeric"), 
                        arrayID = vector(mode = "numeric"), 
                        start = vector(mode = "numeric"), 
                        end = vector(mode = "numeric"), 
                        strand = vector(mode = "character")))
    }

    # export the reference and array sequence
    repeat_file = file.path(cmd_arguments$output_folder, paste0("Array_", i, "_repeat.fasta"))
    array_file = file.path(cmd_arguments$output_folder, paste0("Array_", i, "_sequence.fasta"))
    seqinr::write.fasta(sequences = arrays$representative[i], names = "reference_repeat", file.out = repeat_file, open = "w")
    seqinr::write.fasta(sequences = fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], names = arrays$seqID[i], file.out = array_file, open = "w")
    
    # use nhmmer
    nhmmer_table_output = file.path(cmd_arguments$output_folder, paste0("nhmmer_", i, "_output.txt"))
    system(paste("../dep/hmmer/nhmmer.exe", " ", 
                 "--tblout ", nhmmer_table_output, " ",
                 repeat_file, " ",
                 array_file, sep = ""), 
                 intern = FALSE, wait = TRUE, show.output.on.console = FALSE,  ignore.stderr = TRUE, ignore.stdout = TRUE)
    # read and parse the output into repeat table
    repeats_df <- read_and_format_nhmmer(nhmmer_table_output, arrays$seqID[i], i)
    print(repeats_df)

    if(nrow(repeats_df) == 0) {
      sink()
      gc()
      return(repeats_df)
    }
    # handle overlaps and gaps
    repeats_df = handle_overlaps(repeats_df) # TODO
    print("deb1")
    # adjust start and end coordinates
    repeats_df$start = repeats_df$start + arrays$start[i] - 1
    repeats_df$end = repeats_df$end + arrays$end[i] - 1
    print("deb2")
    # remove generated files
    file.remove(repeat_file)
    file.remove(array_file)
    file.remove(nhmmer_table_output)
    sink()
    gc()
    print("deb3")
    return(repeats_df)
  }
  print("Reporting the output")
  export_gff(annotations.data.frame = repeats,
             output = cmd_arguments$output_folder, 
             file.name = paste0(basename(cmd_arguments$fasta_file), "_repeats"), 
             source = "TRASH",
             type = "Satellite_array", 
             seqid = 1, 
             start = 3, 
             end = 4, 
             attributes = 2, 
             attribute.names = "Array=")

  write.csv(x = repeats, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_repeats.csv")), row.names = FALSE)
  
  stopCluster(cl)
  gc()
}
# Functions TODO:

# split_and_check_arrays: identifies whether a region is a single array or multiple arrays and whether they are simple or complex, returns a list of arrays with these information, maybe N?
#