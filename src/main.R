main <- function(cmd_arguments) {
  cat("TRASH: workspace initialised\n")
  ### Start workers =============================================================================================
  cl <- makeCluster(cmd_arguments$cores_no)
  registerDoParallel(cl)
  foreach (i = 1 : getDoParWorkers()) %dopar% sink() # brings output back to the console

  ### Settings ==================================================================================================
  window_size <- ceiling(cmd_arguments$max_rep_size * 1.1)
  max_eval = 0.001

  ### Load fasta ================================================================================================
  cat(paste0(" Loading the fasta file: ", basename(cmd_arguments$fasta_file), "\n"))
  fasta_content <- read_fasta_and_list(cmd_arguments$fasta_file)
  gc()

  ### Calculate repeat scores for each sequence =================================================================
  cat(" Calculating repeat scores for each sequence\n") # TODO: somewhere on region/array identification step there's a problem with edge telomeric repeats  
  cat("################################################################################\n")
  pb <- txtProgressBar(min = 0, max = length(fasta_content), style = 1)
  repeat_scores <- list()
  for (i in seq_along(fasta_content)) {
    repeat_scores <- append(repeat_scores, list(sequence_window_score(fasta_content[[i]], window_size)))
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
  }
  close(pb)
  gc()

  ### Identify regions with high repeat content and merge into a df =============================================
  cat(" Identifying regions with high repeat content\n")
  repetitive_regions <- data.frame(starts = NULL, ends = NULL, scores = NULL, seqID = NULL, numID = NULL)
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

  ### Split regions into arrays =================================================================================
  cat(" Identifying individual arrays with repeats\n")
  cat("################################################################################\n")
  region_sizes = repetitive_regions$ends - repetitive_regions$starts
  progress_values = region_sizes / sum(region_sizes)
  pb <- txtProgressBar(min = 0, max = 1, style = 1)
  arrays <- foreach (i = seq_len(nrow(repetitive_regions)),
                     .combine = rbind,
                     #.packages = c("foreach", "iterators", "parallel"),
                     .export = c("split_and_check_arrays", "extract_kmers", "collapse_kmers", "genomic_bins_starts", "consensus_N", "write_align_read")) %dopar% {
    #sink()
    out = split_and_check_arrays(start = repetitive_regions$starts[i],
                                  end = repetitive_regions$ends[i],
                                  sequence = fasta_content[[repetitive_regions$numID[i]]][repetitive_regions$starts[i] : repetitive_regions$ends[i]],
                                  seqID = repetitive_regions$seqID[i],
                                  numID = repetitive_regions$numID[i],
                                  arrID = i,
                                  max_repeat = cmd_arguments$max_rep_size,
                                  min_repeat = cmd_arguments$min_rep_size,
                                  mafft = "../dep/mafft-7.520-win64-signed/mafft-win/mafft.bat",
                                  temp_dir = cmd_arguments$output_folder,
                                  src_dir = getwd(),
                                  sink_output = FALSE)
    setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
    return(out)
  }
  close(pb)
  gc()

  ### Save array output =========================================================================================
  cat(" Arrays identified, saving the array table\n") #TODO: move the save to the end of the script and add additional info about the arrays
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
  ### Shift representative repeats and apply templates ==========================================================
  cat(" Shifting array representative repeats and comparing to provided templates\n")
  arrays$representative <- foreach (i = seq_len(nrow(arrays)), .combine = c, .export = "shift_and_compare") %dopar% {
    shift_and_compare(arrays$representative[i], cmd_arguments$templates) #TODO: write it
  }

  ### Map repeats ===============================================================================================
  cat(" Mapping the array representative to the array using nhmmer\n") # TODO: extract this into a single function
  cat("################################################################################\n")
  array_sizes = arrays$end - arrays$start
  progress_values = array_sizes / sum(array_sizes)
  pb <- txtProgressBar(style = 1)
  repeats <- foreach (i = seq_len(nrow(arrays)), .combine = rbind, .export = c("read_and_format_nhmmer", "handle_overlaps", "export_gff")) %dopar% {
    if(arrays$representative[i] == "") {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      return(data.frame(seqID = vector(mode = "numeric"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character")))
    }
    ## export the reference and array sequence ===============================
    repeat_file = file.path(cmd_arguments$output_folder, paste0("Array_", i, "_repeat.fasta"))
    array_file = file.path(cmd_arguments$output_folder, paste0("Array_", i, "_sequence.fasta"))
    seqinr::write.fasta(sequences = arrays$representative[i], names = "reference_repeat", file.out = repeat_file, open = "w")
    seqinr::write.fasta(sequences = fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], names = arrays$seqID[i], file.out = array_file, open = "w")
    ## use nhmmer =============================== ============================
    nhmmer_table_output = file.path(cmd_arguments$output_folder, paste0("nhmmer_", i, "_output.txt"))
    system(paste("../dep/hmmer/nhmmer.exe", " ", 
                 "--dna --tblout ", nhmmer_table_output, " ",
                 repeat_file, " ",
                 array_file, sep = ""), 
                 intern = FALSE, wait = TRUE, show.output.on.console = FALSE, ignore.stderr = FALSE, ignore.stdout = TRUE)
    ## read and parse the output into repeat table ===========================
    repeats_df <- read_and_format_nhmmer(nhmmer_table_output, arrays$seqID[i], i)
    if(nrow(repeats_df) != 0) {
    ## handle overlaps and gaps ==============================================
      repeats_df = handle_overlaps(repeats_df) # TODO, adjust score by -1 for each modification
    ## adjust start and end coordinates ======================================
      repeats_df$start = repeats_df$start + arrays$start[i] - 1
      repeats_df$end = repeats_df$end + arrays$start[i] - 1
    }
    ## remove generated files ================================================
    file.remove(repeat_file)
    file.remove(array_file)
    file.remove(nhmmer_table_output)
    gc()
    setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
    return(repeats_df)
  }
  close(pb)

  ### Write output ==============================================================================================
  cat(" Saving repeats\n")
  repeats$array_N = arrays$top_N[repeats$arrayID]
  repeats$score_by_len = repeats$score / repeats$array_N
  repeats = repeats[repeats$eval <= max_eval,]
  export_gff(annotations.data.frame = repeats,
             output = cmd_arguments$output_folder, 
             file.name = paste0(basename(cmd_arguments$fasta_file), "_repeats"), 
             source = "TRASH",
             type = "Satellite_DNA", 
             seqid = 1, 
             start = 3, 
             end = 4, 
             strand = 5,
             attributes = c(2,9,7), 
             attribute.names = c("Array=", "Score=", "Eval="))

  write.csv(x = repeats, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_repeats.csv")), row.names = FALSE)
  #str(repeats)
  #str(arrays)
  stopCluster(cl)
  gc()
}
