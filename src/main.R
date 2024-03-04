main <- function(cmd_arguments) {
  cat("TRASH: workspace initialised                                     ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  # TODO: remove this development settings
  if(Sys.info()['sysname'] == "Windows") {
    mafft_dir = "../dep/mafft-7.520-win64-signed/mafft-win/mafft.bat"
    # nhmmer_dir = "../dep/hmmer/nhmmer.exe"
    nhmmer_dir = "C:/cygwin64/home/Piotr Włodzimierz/hmmer/hmmer-3.4/src/nhmmer.exe"
  } else {
    mafft_dir = "mafft"
    nhmmer_dir = "nhmmer"
  }

  ### 01 / 14 Start workers =============================================================================================
  log_messages <- file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_main_log_file.txt")) #TODO make into a flag

  cl <- makeCluster(cmd_arguments$cores_no)
  registerDoParallel(cl)
  foreach (i = 1 : getDoParWorkers()) %dopar% sink() # brings output back to the console

  if (log_messages != "") {
    cat("Log messages of file:  ", basename(cmd_arguments$fasta_file), "\n",
        file = log_messages, append = FALSE)
    cat("\nTime:        ", date(), "\n",
        file = log_messages, append = TRUE)
    cat("Max repeat size: ", cmd_arguments$max_rep_size, "\nMin repeat size: ", cmd_arguments$min_rep_size, "\nTemplates file: ", cmd_arguments$templates, "\nCores number: ", cmd_arguments$cores_no, "\n",
        file = log_messages, append = TRUE)
  }

  ### 02 / 14 Settings ==================================================================================================
  window_size <- ceiling(cmd_arguments$max_rep_size * 1.1)
  use_adist_scores = TRUE # using nhmmer, recalculate scores for consistency with other methods
  fix_overlaps = TRUE
  fix_gaps = TRUE
  do_shift_classes <- TRUE

  ### 03 / 14 Load fasta ================================================================================================
  cat(paste0(" 03 / 13 Loading the fasta file: ", basename(cmd_arguments$fasta_file)))
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  fasta_content <- read_fasta_and_list(cmd_arguments$fasta_file)
  gc()
  if (log_messages != "") cat("03 / 14 \nFasta sizes: ", sapply(fasta_content, length), "\n", file = log_messages, append = TRUE)

  ### 04 / 14 Calculate repeat scores for each sequence =================================================================
  cat(" 04 / 13 Calculating repeat scores for each sequence             ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  pb <- txtProgressBar(min = 0, max = length(fasta_content), style = 1)
  repeat_scores <- list()
  for (i in seq_along(fasta_content)) {
    repeat_scores <- append(repeat_scores, list(sequence_window_score(fasta_content[[i]], window_size, log_messages))) # it's parallel
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
  }
  close(pb)
  gc()
  if (log_messages != "") cat("04 / 14 \nFinished scores, length of the list: ", sapply(repeat_scores, length), "\n", file = log_messages, append = TRUE)

  ### 05 / 14 Identify regions with high repeat content and merge into a df =============================================
  cat(" 05 / 13 Identifying regions with high repeat content            ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  repetitive_regions <- data.frame(starts = NULL, ends = NULL, scores = NULL, seqID = NULL, numID = NULL)
  for (i in seq_along(repeat_scores)) {
    # if (log_messages != "") cat(" ", i, "\n", file = log_messages, append = TRUE)
    # if (log_messages != "") cat(" ", length(fasta_content[[i]]), "\n", file = log_messages, append = TRUE)
    # if (log_messages != "") cat(" ", length(repeat_scores[[i]]), "\n", file = log_messages, append = TRUE)
    regions_of_sequence <- merge_windows(list_of_scores = repeat_scores[[i]], window_size = window_size, sequence_full_length = length(fasta_content[[i]]), log_messages)
    if (nrow(regions_of_sequence) != 0) {
      regions_of_sequence$seqID <- names(fasta_content)[[i]]
      regions_of_sequence$numID <- i
      repetitive_regions <- rbind(repetitive_regions, regions_of_sequence)
    }
  }
  if(!inherits(repetitive_regions, "data.frame")) stop("No regions with repeats identified")
  if(nrow(repetitive_regions) == 0) stop("No regions with repeats identified")

  if (log_messages != "") cat("\n05 / 14 \nFinished regions, nrow of the regions: ", nrow(repetitive_regions), "\n", file = log_messages, append = TRUE)

  write.csv(x = repetitive_regions, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrays.csv")), row.names = FALSE) # TODO remove

  ### 06 / 14 Split regions into arrays =================================================================================
  cat(" 06 / 13 Identifying individual arrays with repeats              ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  region_sizes = repetitive_regions$ends - repetitive_regions$starts
  progress_values = region_sizes / sum(region_sizes)
  pb <- txtProgressBar(min = 0, max = 1, style = 1)
  if (log_messages != "") cat("\n06 / 14 Where's this error?", file = log_messages, append = TRUE)
  arrays <- foreach (i = seq_len(nrow(repetitive_regions)),
                     .combine = rbind,
                     .export = c("split_and_check_arrays", "extract_kmers", "collapse_kmers", "genomic_bins_starts", "consensus_N", "write_align_read")) %dopar% {
    # if (log_messages != "") cat("\n Got in region", i, file = log_messages, append = TRUE)
    #.libPaths(c(.libPaths(), gsub("src", "R_libs", getwd())))
    out = split_and_check_arrays(start = repetitive_regions$starts[i],
                                  end = repetitive_regions$ends[i],
                                  sequence = fasta_content[[repetitive_regions$numID[i]]][repetitive_regions$starts[i] : repetitive_regions$ends[i]],
                                  seqID = repetitive_regions$seqID[i],
                                  numID = repetitive_regions$numID[i],
                                  arrID = i,
                                  max_repeat = cmd_arguments$max_rep_size,
                                  min_repeat = cmd_arguments$min_rep_size,
                                  mafft = mafft_dir,
                                  temp_dir = cmd_arguments$output_folder,
                                  src_dir = getwd(),
                                  sink_output = FALSE) # TODO: Extract N calc from array splitting
    setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
    return(out)
  }
  close(pb)
  gc()

  if (log_messages != "") cat("06 / 14 \nFinished arrays, nrow of the arrays: ", nrow(arrays), "\n", file = log_messages, append = TRUE)

  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrays.csv")), row.names = FALSE) # TODO remove

  ### 07 / 14 Shift representative repeats and apply templates ==========================================================
  cat(" 07 / 13 Shifting representative and comparing templates         ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  pb <- txtProgressBar(min = 0, max = nrow(arrays), style = 1)
  if (log_messages != "") cat("07; set templates\n", file = log_messages, append = TRUE)
  if(cmd_arguments$templates != 0) {
    templates = read_fasta_and_list(cmd_arguments$templates)
    if(length(templates) == 0) stop("No templates found within the template file")
  } else {
    templates = 0
  }
  if (log_messages != "") cat("07; start representative shift\n", file = log_messages, append = TRUE)
  arrays$representative <- foreach (i = seq_len(nrow(arrays)), .combine = c, .export = c("shift_and_compare", "shift_sequence", "compare_circular", "rev_comp_string", "kmer_hash_score")) %dopar% {
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    # if (log_messages != "") cat("07; shifting ", arrays$representative[i],"\n", file = log_messages, append = TRUE)
    if(!inherits(arrays$representative[i], "character")) return("_")
    return(shift_and_compare(arrays$representative[i], templates))
  }
  if (log_messages != "") cat("07; finished representative shift\n", file = log_messages, append = TRUE)
  arrays$class = ""
  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrdays.csv")), row.names = FALSE) # TODO remove

  if (log_messages != "") cat("07; assigning classes\n", file = log_messages, append = TRUE)
  for(i in seq_len(nrow(arrays))) {
    if(arrays$representative[i] == "_") {arrays$representative[i] = ""; next()}
    arrays$class[i] = strsplit(arrays$representative[i], split = "_")[[1]][1]
    arrays$representative[i] = strsplit(arrays$representative[i], split = "_")[[1]][2]
  }
  close(pb)

  if (log_messages != "") cat("07 / 14 \nFinished representative, nrow of the arrays: ", nrow(arrays), "\n", file = log_messages, append = TRUE)

  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrays.csv")), row.names = FALSE) # TODO remove

  ### 08 / 14 Classify unclassified and shift ===========================================================================
  cat(" 08 / 13 Classifying remaining representative repeats            ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  if (log_messages != "") cat("8 class start: ", Sys.time(), "\n", file = log_messages, append = TRUE)
  arrays = classify_repeats(repeat_df = arrays)
  ## Shift representatives to match the "most important" one ================
  if(do_shift_classes) {
    classes <- unique(arrays$class)
    classes <- classes[!(classes %in% c(names(templates), "none_identified"))]
    if (log_messages != "") cat("8 shift start: ", Sys.time(), "\n", file = log_messages, append = TRUE)
    arrays_t <- foreach (i = seq_along(classes), .combine = rbind, .export = c("shift_classes", "compare_circular", "rev_comp_string")) %dopar% {
      # if (log_messages != "") cat("08 / 14 \nClassifying class ", i, "\n", file = log_messages, append = TRUE)
      arrays_class <- arrays[arrays$class == classes[i], ]
      arrays_class$representative <- shift_classes(arrays_class)
      return(arrays_class)
    }
    arrays <- rbind(arrays_t, arrays[which(arrays$class %in% c(names(templates), "none_identified")), ])
    remove(arrays_t)
  }
  if (log_messages != "") cat("8 finito: ", Sys.time(), "\n", file = log_messages, append = TRUE)

  if (log_messages != "") cat("08 / 14 \nFinished classifications, nrow of the arrays: ", nrow(arrays), "\n", file = log_messages, append = TRUE)

  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrays.csv")), row.names = FALSE) # TODO remove
  ### 09 / 14 Map repeats ===============================================================================================
  cat(" 09 / 13 Mapping the array representative to the array           ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  array_sizes = arrays$end - arrays$start
  progress_values = array_sizes / sum(array_sizes)
  pb <- txtProgressBar(style = 1)
  if (log_messages != "") cat("09; map loop\n", file = log_messages, append = TRUE)
  repeats <- foreach (i = seq_len(nrow(arrays)), .combine = rbind, .export = c("write_align_read", "consensus_N", "read_and_format_nhmmer", "handle_overlaps", "handle_gaps", "export_gff", "map_nhmmer", "map_default", "rev_comp_string")) %dopar% {
    time_report_df = as.numeric(Sys.time())
    if (arrays$representative[i] == "") {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      return(data.frame(seqID = vector(mode = "character"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric"),
                        eval = vector(mode = "numeric"),
                        width = vector(mode = "numeric"),
                        class = vector(mode = "character"),
                        representative = vector(mode = "character"),
                        score_template = vector(mode = "numeric")))
    }
    if(arrays$top_N[i] >= 14) {
    ## nhmmer for repeats of 14+ bp =========================================
      repeats_df = map_nhmmer(cmd_arguments$output_folder, i, arrays$representative[i], arrays$seqID[i], arrays$start[i], 
                              arrays$end[i], fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], nhmmer_dir)
    } else {
    ## matchpattern for shorter =============================================
      repeats_df = map_default(i, arrays$representative[i], arrays$seqID[i], arrays$start[i], paste(fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], collapse = ""))
    }
    time_report_df = c(time_report_df, as.numeric(Sys.time()))
    
    if(nrow(repeats_df) < 2) {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      return(data.frame(seqID = vector(mode = "character"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric"),
                        eval = vector(mode = "numeric"),
                        width = vector(mode = "numeric"),
                        class = vector(mode = "character"),
                        representative = vector(mode = "character"),
                        score_template = vector(mode = "numeric")))
    }
    ## add width and class ============================================================
    repeats_df$width = repeats_df$end - repeats_df$start + 1
    repeats_df$class = arrays$class[i]
    ## handle overlaps ======================================================
    
    if (fix_overlaps) repeats_df = handle_overlaps(repeats_df, overlap_threshold = 0.1, representative_len = arrays$top_N[i])
    if(nrow(repeats_df) < 2) {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      return(data.frame(seqID = vector(mode = "character"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric"),
                        eval = vector(mode = "numeric"),
                        width = vector(mode = "numeric"),
                        class = vector(mode = "character"),
                        representative = vector(mode = "character"),
                        score_template = vector(mode = "numeric")))
    }
    time_report_df = c(time_report_df, as.numeric(Sys.time()))
    ## handle gaps if proper array ==========================================
    
    if (fix_gaps) {
      repeats_df <- handle_gaps(repeats_df, representative_len = arrays$top_N[i])
      # if (sum(repeats_df$width) > (((arrays$end[i] - arrays$start[i]) - cmd_arguments$max_rep_size * 2) / 2)) {
      #   repeats_df <- handle_gaps(repeats_df, representative_len = arrays$top_N[i])
      # } else {
      #   for (j in seq_len(nrow(repeats_df))) repeats_df$class[j] <- paste0(repeats_df$class[j], "_scattered")
      # }
    }
    time_report_df = c(time_report_df, as.numeric(Sys.time()))
    
    ## Return nothing if handle_gaps removed all repeats ====================
    if (nrow(repeats_df) < 2) {
      # if (log_messages != "" && i == 1) cat("09; array no ", i, " no repeats after gaps\n", file = log_messages, append = TRUE)
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      return(data.frame(seqID = vector(mode = "character"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric"),
                        eval = vector(mode = "numeric"),
                        width = vector(mode = "numeric"),
                        class = vector(mode = "character"),
                        representative = vector(mode = "character"),
                        score_template = vector(mode = "numeric")))
    }
    ## Recalculate representative ===========================================
    # TODO: Use more repeats for long arrays, this is stringent and good for most arrays, but some deserve a better recalculation
    # TODO: check if its the mafft call that is the long part here and see how write_align_read can be sped up
    max_repeats_to_align <- 15
    min_repeats_to_recalculate <- 10
    repeats_df$representative <- arrays$representative[i]
    repeats_df$score_template <- -1
    sample_IDs = which(repeats_df$strand != ".")
    if (length(sample_IDs) >= min_repeats_to_recalculate) {
      if (length(sample_IDs) > max_repeats_to_align) sample_IDs <- sample(sample_IDs, max_repeats_to_align)
      repeats_seq = unlist(lapply(sample_IDs, function(X) paste0(fasta_content[[arrays$numID[i]]][repeats_df$start[X] : repeats_df$end[X]], collapse = "")))
      strands = repeats_df$strand[sample_IDs]
      repeats_seq[which(strands == "-")] = unlist(lapply(repeats_seq[which(strands == "-")], rev_comp_string))
      alignment <- write_align_read(mafft_exe = mafft_dir,
                                    temp_dir = cmd_arguments$output_folder,
                                    sequences = repeats_seq,
                                    name = paste(basename(cmd_arguments$fasta_file), arrays$seqID[i], arrays$numID[i], i, runif(1, 0, 1), sep = "_"))
      consensus <- consensus_N(alignment, arrays$top_N[i])
      if (length(consensus) != 0) repeats_df$representative <- consensus
    }
    time_report_df = c(time_report_df, as.numeric(Sys.time()))

    # if (log_messages != "") cat("09; array no ", i, " representative finished\n", file = log_messages, append = TRUE)

    ## If set, change to edit distance based score =========================
    # TODO: use only unique repeats to recalculate, should make it run faster for arrays with many repeats (where many are also identical)
    if(use_adist_scores) {
      repeats_seq = unlist(lapply(seq_len(nrow(repeats_df)), function(X) paste0(fasta_content[[arrays$numID[i]]][repeats_df$start[X] : repeats_df$end[X]], collapse = "")))
      costs = list(insertions = 1, deletions = 1, substitutions = 1)
      if (sum(repeats_df$strand == "+") > 0) repeats_df$score[repeats_df$strand == "+"] = adist(repeats_df$representative[1], repeats_seq[repeats_df$strand == "+"], costs)[1,]  / nchar(repeats_df$representative[1]) * 100
      if (sum(repeats_df$strand == "-") > 0) repeats_df$score[repeats_df$strand == "-"] = adist(rev_comp_string(repeats_df$representative[1]), repeats_seq[repeats_df$strand == "-"])[1,]  / nchar(repeats_df$representative[1]) * 100

      if (arrays$class[i] %in% names(templates)) {
        template <- paste(templates[[which(arrays$class[i] == names(templates))]], collapse = "")
        if (sum(repeats_df$strand == "+") > 0) repeats_df$score_template[repeats_df$strand == "+"] = adist(template, repeats_seq[repeats_df$strand == "+"], costs)[1,]  / nchar(template) * 100
        if (sum(repeats_df$strand == "-") > 0) repeats_df$score_template[repeats_df$strand == "-"] = adist(rev_comp_string(template), repeats_seq[repeats_df$strand == "-"])[1,]  / nchar(template) * 100
      }
    }
    time_report_df = c(time_report_df, as.numeric(Sys.time()))
    gc()
    setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
    if (nrow(repeats_df) < 2) {
      return(data.frame(seqID = vector(mode = "character"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric"),
                        eval = vector(mode = "numeric"),
                        width = vector(mode = "numeric"),
                        class = vector(mode = "character"),
                        representative = vector(mode = "character"),
                        score_template = vector(mode = "numeric")))
    }
    repeats_df <- repeats_df[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class", "representative", "score_template")]
    repeats_df$seqID = as.character(repeats_df$seqID)
    repeats_df$arrayID = as.numeric(repeats_df$arrayID)
    repeats_df$start = as.numeric(repeats_df$start)
    repeats_df$end = as.numeric(repeats_df$end)
    repeats_df$strand = as.character(repeats_df$strand)
    repeats_df$score = as.numeric(repeats_df$score)
    repeats_df$eval = as.numeric(repeats_df$eval)
    repeats_df$width = as.numeric(repeats_df$width)
    repeats_df$class = as.character(repeats_df$class)
    repeats_df$representative = as.character(repeats_df$representative)
    repeats_df$score_template = as.numeric(repeats_df$score_template)
    return(repeats_df)
  }
  close(pb)

  if (log_messages != "") cat("09 / 14 \nFinished mapping, nrow of the repeats: ", nrow(repeats), "\n", file = log_messages, append = TRUE)

  ### 10 / 14 ===========================================================================================================
  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrays.csv")), row.names = FALSE)

  ### 11 / 14 Summarise array information ===============================================================================
  cat(" 11 / 13 Summarising array information                           ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")


  # # Find unique array IDs in repeats
  # unique_ids <- unique(repeats$arrayID)

  # # Subset repeats to only include unique IDs
  # unique_repeats <- repeats[match(unique_ids, repeats$arrayID), ]

  # # Match array IDs in arrays with unique IDs in repeats
  # match_indices <- match(arrays$arrayID, unique_ids)

  # # Identify non-NA indices
  # non_na_indices <- !is.na(match_indices)

  # # Update representative values in arrays where match is found
  # arrays$representative[non_na_indices] <- unique_repeats$representative[match_indices[non_na_indices]]

  for(i in seq_len((nrow(arrays)))) {
    if(sum((repeats$arrayID == i) > 0)) {
      arrays$representative[i] <- repeats$representative[which(repeats$arrayID == i)[1]]
    }
  }
  repeats <- repeats[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class", "score_template")]


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

  if (log_messages != "") cat("10 / 14 \nFinished summarising, nrow of the arrays: ", nrow(arrays), "\n", file = log_messages, append = TRUE)

  ### 12 / 14 Save array output =========================================================================================
  cat(" 12 / 13 Saving the array table                                  ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
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
  cat(" 13 / 13 Saving the repeats table                                ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
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

  if (log_messages != "") cat("14 / 14 \n## Done ##\nTime:         ", date(), "\n", file = log_messages, append = TRUE)

  stopCluster(cl)
  gc()
}
