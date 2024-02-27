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
  cl <- makeCluster(cmd_arguments$cores_no)
  registerDoParallel(cl)
  foreach (i = 1 : getDoParWorkers()) %dopar% sink() # brings output back to the console

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

  ### 04 / 14 Calculate repeat scores for each sequence =================================================================
  cat(" 04 / 13 Calculating repeat scores for each sequence             ")
  cat(Sys.time())
  cat("\n")
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
  cat(" 05 / 13 Identifying regions with high repeat content            ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
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
  cat(" 06 / 13 Identifying individual arrays with repeats              ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  region_sizes = repetitive_regions$ends - repetitive_regions$starts
  progress_values = region_sizes / sum(region_sizes)
  pb <- txtProgressBar(min = 0, max = 1, style = 1)
  arrays <- foreach (i = seq_len(nrow(repetitive_regions)),
                     .combine = rbind,
                     .export = c("split_and_check_arrays", "extract_kmers", "collapse_kmers", "genomic_bins_starts", "consensus_N", "write_align_read")) %dopar% {
    .libPaths(c(.libPaths(), gsub("src", "R_libs", getwd())))
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

  ### 07 / 14 Shift representative repeats and apply templates ==========================================================
  cat(" 07 / 13 Shifting representative and comparing templates         ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  pb <- txtProgressBar(min = 0, max = nrow(arrays), style = 1)
  if(cmd_arguments$templates != "") {
    templates = read_fasta_and_list(templates)
    if(length(templates == 0)) stop("No templates found within the template file")
  } else templates = ""
  arrays$representative <- foreach (i = seq_len(nrow(arrays)), .combine = c, .export = c("shift_and_compare", "shift_sequence", "compare_circular", "rev_comp_string", "kmer_hash_score")) %dopar% {
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    return(shift_and_compare(arrays$representative[i], templates))
  }
  arrays$class = ""
  for(i in seq_len(nrow(arrays))) {
    if(arrays$representative[i] == "_") {arrays$representative[i] = ""; next()}
    arrays$class[i] = strsplit(arrays$representative[i], split = "_")[[1]][1]
    arrays$representative[i] = strsplit(arrays$representative[i], split = "_")[[1]][2]
  }
  close(pb)

  ### 08 / 14 Classify unclassified and shift ===========================================================================
  cat(" 08 / 13 Classifying remaining representative repeats            ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  arrays = classify_repeats(repeat_df = arrays)
  ## Shift representatives to match the "most important" one ================
  if(do_shift_classes) {
    classes <- unique(arrays$class)
    classes <- classes[!(classes %in% c(names(templates), "none_identified"))]
    arrays_t <- foreach (i = seq_along(classes), .combine = rbind, .export = c("shift_classes", "compare_circular", "rev_comp_string")) %dopar% {
      arrays_class <- arrays[arrays$class == classes[i], ]
      arrays_class$representative <- shift_classes(arrays_class)
      return(arrays_class)
    }
    arrays <- rbind(arrays_t, arrays[!(arrays$classes %in% classes)])
    remove(arrays_t)
  }


  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrays.csv")), row.names = FALSE) # TODO remove
  ### 09 / 14 Map repeats ===============================================================================================
  cat(" 09 / 13 Mapping the array representative to the array           ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  array_sizes = arrays$end - arrays$start
  progress_values = array_sizes / sum(array_sizes)
  pb <- txtProgressBar(style = 1)
  repeats <- foreach (i = seq_len(nrow(arrays)), .combine = rbind, .export = c("write_align_read", "consensus_N", "read_and_format_nhmmer", "handle_overlaps", "handle_gaps", "export_gff", "map_nhmmer", "map_default", "rev_comp_string")) %dopar% {
    #cat(paste0("Started ", i, "\n"), file = file.path(cmd_arguments$output_folder, paste0(i,"log_temp.txt")), append = TRUE)
    if(arrays$representative[i] == "") {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      #cat(paste0("Finished ", i, "\n"), file = file.path(cmd_arguments$output_folder, paste0(i,"log_temp.txt")), append = TRUE)
      return(data.frame(seqID = vector(mode = "numeric"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric",),
                        eval = vector(mode = "numeric"),
                        class = vector(mode = "character")))
    }
    if(arrays $top_N[i] >= 14) {
    ## nhmmer for repeats of 14+ bp =========================================
      repeats_df = map_nhmmer(cmd_arguments$output_folder, i, arrays$representative[i], arrays$seqID[i], arrays$start[i], 
                              arrays$end[i], fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], nhmmer_dir)
    } else {
    ## matchpattern for shorter =============================================
      repeats_df = map_default(i, arrays$representative[i], arrays$seqID[i], arrays$start[i], paste(fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], collapse = ""))
    }
    if(nrow(repeats_df) == 0) {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      #cat(paste0("Finished ", i, "\n"), file = file.path(cmd_arguments$output_folder, paste0(i,"log_temp.txt")), append = TRUE)
      return(data.frame(seqID = vector(mode = "numeric"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric",),
                        eval = vector(mode = "numeric"),
                        class = vector(mode = "character")))
    }
    #cat(paste0(nrow(repeats_df), "\n"), file = file.path(cmd_arguments$output_folder, paste0(i,"log_temp.txt")), append = TRUE)
    ## add width ============================================================
    repeats_df$width = repeats_df$end - repeats_df$start + 1
    ## handle overlaps ======================================================
    if (fix_overlaps && (nrow(repeats_df) > 1)) repeats_df = handle_overlaps(repeats_df, overlap_threshold = 0.1, representative_len = arrays$top_N[i])
    ## handle gaps if proper array ==========================================
    repeats_df$class = arrays$class[i]
    if (fix_gaps && (nrow(repeats_df) > 1)) {
      if (sum(repeats_df$width) > (((arrays$end[i] - arrays$start[i]) - cmd_arguments$max_rep_size * 2) / 2)) {
        repeats_df <- handle_gaps(repeats_df, overlap_threshold = 0.1, representative_len = arrays$top_N[i])
      } else {
        for (j in seq_len(nrow(repeats_df))) repeats_df$class[j] <- paste0(repeats_df$class[j], "_scattered")
      }
    }
    #cat(paste0(nrow(repeats_df), "\n"), file = file.path(cmd_arguments$output_folder, paste0(i,"log_temp.txt")), append = TRUE)
    ## Return nothing if handle_gaps removed all repeats ====================
    if (nrow(repeats_df) == 0) {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      #cat(paste0("Finished ", i, "\n"), file = file.path(cmd_arguments$output_folder, paste0(i,"log_temp.txt")), append = TRUE)
      return(data.frame(seqID = vector(mode = "numeric"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric",),
                        eval = vector(mode = "numeric"),
                        class = vector(mode = "character")))
    }
    ## Recalculate representative ===========================================
    max_repeats_to_align <- 50
    min_repeats_to_recalculate <- 5
    sample_IDs = which(repeats_df$strand != ".")
    if(length(sample_IDs) >= min_repeats_to_recalculate) {
      if(length(sample_IDs) > max_repeats_to_align) sample_IDs <- sample(sample_IDs, max_repeats_to_align)
      repeats_seq = unlist(lapply(sample_IDs, function(X) paste0(fasta_content[[arrays$numID[i]]][repeats_df$start[X] : repeats_df$end[X]], collapse = "")[[1]]))
      strands = repeats_df$strand[sample_IDs]
      repeats_seq[which(strands == "-")] = unlist(lapply(repeats_seq[which(strands == "-")], rev_comp_string))
      alignment <- write_align_read(mafft_exe = mafft_dir,
                                    temp_dir = cmd_arguments$output_folder,
                                    sequences = repeats_seq,
                                    name = paste(arrays$seqID[i], arrays$numID[i], i, runif(1, 0, 1), sep = "_"))
      consensus <- consensus_N(alignment, arrays$top_N[i])
      if(length(consensus) != 0) arrays$representative[i] <- consensus
    }
    #cat(paste0(nrow(repeats_df), "\n"), file = file.path(cmd_arguments$output_folder, paste0(i,"log_temp.txt")), append = TRUE)

    ## If set, change to edit distance based score =========================
    if(use_adist_scores) {
      repeats_seq = unlist(lapply(seq_len(nrow(repeats_df)), function(X) paste0(fasta_content[[arrays$numID[i]]][repeats_df$start[X] : repeats_df$end[X]], collapse = "")[[1]]))
      costs = list(insertions = 1, deletions = 1, substitutions = 1)
      repeats_df$score[repeats_df$strand == "+"] = adist(arrays$representative[i], repeats_seq[repeats_df$strand == "+"], costs)[1,]  / nchar(arrays$representative[i]) * 100
      repeats_df$score[repeats_df$strand == "-"] = adist(rev_comp_string(arrays$representative[i]), repeats_seq[repeats_df$strand == "-"])[1,]  / nchar(arrays$representative[i]) * 100
    }
    #cat(paste0(nrow(repeats_df), "\n"), file = file.path(cmd_arguments$output_folder, paste0(i,"log_temp.txt")), append = TRUE)

    gc()
    setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
    #cat(paste0("Finished ", i, "\n"), file = file.path(cmd_arguments$output_folder, paste0(i,"log_temp.txt")), append = TRUE)
    return(repeats_df)
  }
  close(pb)
  ### 10 / 14 ===========================================================================================================
  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrays.csv")), row.names = FALSE)

  ### 11 / 14 Summarise array information ===============================================================================
  cat(" 11 / 13 Summarising array information                           ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
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
  stopCluster(cl)
  gc()
}
