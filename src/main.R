main <- function(cmd_arguments) {
  cat("TRASH: workspace initialised\n")

  ### 01 / 14 Start workers =============================================================================================
  cl <- makeCluster(cmd_arguments$cores_no)
  registerDoParallel(cl)
  foreach (i = 1 : getDoParWorkers()) %dopar% sink() # brings output back to the console

  ### 02 / 14 Settings ==================================================================================================
  window_size <- ceiling(cmd_arguments$max_rep_size * 1.1)
  max_eval = 0.001
  use_adist_scores = TRUE # using nhmmer, recalculate scores for consistency with other methods
  fix_overlaps = TRUE
  fix_gaps = TRUE

  ### 03 / 14 Load fasta ================================================================================================
  cat(paste0(" Loading the fasta file: ", basename(cmd_arguments$fasta_file), "\n"))
  fasta_content <- read_fasta_and_list(cmd_arguments$fasta_file)
  gc()

  ### 04 / 14 Calculate repeat scores for each sequence =================================================================
  cat(" Calculating repeat scores for each sequence\n")
  cat("################################################################################\n")
  pb <- txtProgressBar(min = 0, max = length(fasta_content), style = 1)
  repeat_scores <- list()
  for (i in seq_along(fasta_content)) {
    repeat_scores <- append(repeat_scores, list(sequence_window_score(fasta_content[[i]], window_size)))
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
  }
  close(pb)
  gc()

  ### 05 / 14 Identify regions with high repeat content and merge into a df =============================================
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

  ### 06 / 14 Split regions into arrays =================================================================================
  cat(" Identifying individual arrays with repeats\n")
  cat("################################################################################\n")
  region_sizes = repetitive_regions$ends - repetitive_regions$starts
  progress_values = region_sizes / sum(region_sizes)
  pb <- txtProgressBar(min = 0, max = 1, style = 1)
  arrays <- foreach (i = seq_len(nrow(repetitive_regions)),
                     .combine = rbind,
                     .export = c("split_and_check_arrays", "extract_kmers", "collapse_kmers", "genomic_bins_starts", "consensus_N", "write_align_read")) %dopar% {
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
                                  sink_output = TRUE)
    setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
    return(out)
  }
  close(pb)
  gc()

  ### 07 / 14 Shift representative repeats and apply templates ==========================================================
  cat(" Shifting array representative repeats and comparing to provided templates\n")
  arrays$representative <- foreach (i = seq_len(nrow(arrays)), .combine = c, .export = "shift_and_compare") %dopar% {
    shift_and_compare(arrays$representative[i], cmd_arguments$templates) #TODO: make it
  }

  ### 08 / 14 Map repeats ===============================================================================================
  cat(" Mapping the array representative to the array\n")
  cat("################################################################################\n")
  array_sizes = arrays$end - arrays$start
  progress_values = array_sizes / sum(array_sizes)
  pb <- txtProgressBar(style = 1)
  repeats <- foreach (i = seq_len(nrow(arrays)), .combine = rbind, .export = c("read_and_format_nhmmer", "handle_overlaps", "handle_gaps", "export_gff", "map_nhmmer", "map_default", "rev_comp_string")) %dopar% {
    if(arrays$representative[i] == "") {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      return(data.frame(seqID = vector(mode = "numeric"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric",),
                        eval = vector(mode = "numeric")))
    }
    if(arrays$top_N[i] >= 14) {
    ## nhmmer for repeats of 14+ bp =========================================
      repeats_df = map_nhmmer(cmd_arguments$output_folder, i, arrays$representative[i], arrays$seqID[i], arrays$start[i], 
                              arrays$end[i], fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], use_adist_scores)
    } else {
    ## matchpattern for shorter =============================================
      repeats_df = map_default(i, arrays$representative[i], arrays$seqID[i], arrays$start[i], paste(fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], collapse = ""))
    }
    ## add width ============================================================
    repeats_df$width = repeats_df$end - repeats_df$start + 1
    ## handle overlaps and gaps =============================================
    if(fix_overlaps && (nrow(repeats_df) > 1)) repeats_df = handle_overlaps(repeats_df, overlap_threshold = 0.1, representative_len = arrays$top_N[i])
    if(fix_gaps && (nrow(repeats_df) > 1)) repeats_df = handle_gaps(repeats_df, overlap_threshold = 0.1, representative_len = arrays$top_N[i])
  
    gc()
    setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
    return(repeats_df)
  }
  close(pb)

  ### 09 / 14 Filter out low eval score repeats =========================================================================
  # repeats = repeats[repeats$eval <= max_eval,] # not good for short repeats

  ### 10 / 14 ===========================================================================================================

  ### 11 / 14 Summarise array information ===============================================================================
  arrays$repeats_number = 0
  arrays$median_repeat_width = 0
  arrays$median_score = -1

  for (i in seq_len(nrow(arrays))) {
    repeats_temp = repeats[repeats$arrayID == i,]
    if(nrow(repeats_temp) > 0) {
      arrays$repeats_number[i] <- nrow(repeats_temp)
      arrays$median_repeat_width[i] <- ceiling(median(repeats_temp$width))
      if(use_adist_scores) arrays$median_score[i] <- ceiling(median(repeats_temp$score))
    }
  }

  ### 12 / 14 Save array output =========================================================================================
  cat(" Saving the array table\n") #TODO: move the save to the end of the script and add additional info about the arrays
  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrays.csv")), row.names = FALSE)
  export_gff(annotations.data.frame = arrays,
             output = cmd_arguments$output_folder, 
             file.name = paste0(basename(cmd_arguments$fasta_file), "_arrays"), 
             source = "TRASH",
             type = "Satellite_array", 
             seqid = 3, 
             start = 1, 
             end = 2, 
             score = 5,
             attributes = 6, 
             attribute.names = "Name=")

  ### 13 / 14 Save repeat output ========================================================================================
  cat(" Saving the repeats table\n")
  write.csv(x = repeats, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_repeats.csv")), row.names = FALSE)
  export_gff(annotations.data.frame = repeats,
             output = cmd_arguments$output_folder, 
             file.name = paste0(basename(cmd_arguments$fasta_file), "_repeats"), 
             source = "TRASH",
             type = "Satellite_DNA", 
             seqid = 1, 
             start = 3, 
             end = 4, 
             strand = 5,
             attributes = c(2,6,7), 
             attribute.names = c("Array=", "Score=", "Eval="))
  
  ### 14 / 14 Done ======================================================================================================
  stopCluster(cl)
  gc()
}
