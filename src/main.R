main <- function(cmd_arguments) {
  cat("TRASH: workspace initialised                                     ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  # TODO: remove this development settings
  if (Sys.info()["sysname"] == "Windows") {
    mafft_dir <- "../dep/mafft-7.520-win64-signed/mafft-win/mafft.bat"
    # nhmmer_dir <- "../dep/hmmer/nhmmer.exe"
    nhmmer_dir <- "C:/cygwin64/home/Piotr Włodzimierz/hmmer/hmmer-3.4/src/nhmmer.exe"
  } else {
    mafft_dir <- "mafft"
    nhmmer_dir <- "nhmmer"
  }
  cat("================================================================================\n")

  ### 01 / 14 Start workers =============================================================================================
  log_messages <- file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_main_log_file.txt")) #TODO make into a flag
  log_messages <- ""

  cl <- makeCluster(cmd_arguments$cores_no)
  registerDoParallel(cl)
  foreach (i = 1 : getDoParWorkers()) %dopar% {
    sink() # brings output back to the console
    set.seed(0) # Sets random seed for reproducibility
  }

  if (log_messages != "") {
    cat("Log messages of file:  ", basename(cmd_arguments$fasta_file), "\n",
        file = log_messages, append = FALSE)
    cat("\nTime:        ", date(), "\n",
        file = log_messages, append = TRUE)
    cat("Max repeat size: ", cmd_arguments$max_rep_size, "\nMin repeat size: ", cmd_arguments$min_rep_size, "\nTemplates file: ", cmd_arguments$templates, "\nCores number: ", cmd_arguments$cores_no, "\n",
        file = log_messages, append = TRUE)
  }

  ### 02 / 14 Settings ==================================================================================================
  kmer <- 10
  window_size <- (cmd_arguments$max_rep_size + kmer) * 2
  report_runtime <- TRUE

  times <- list(time = as.numeric(Sys.time()), event = "Start main function", data_type = "none", data_value = 0)

  ### 03 / 14 Load fasta ================================================================================================
  if (nchar(basename(cmd_arguments$fasta_file)) > 32) spaces <- 1 else spaces <- 32 - nchar(basename(cmd_arguments$fasta_file))
  cat(paste0(" 03 / 13 Loading the fasta file: ", basename(cmd_arguments$fasta_file), paste(rep(" ", spaces), collapse = "")))
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  fasta_content <- read_fasta_and_list(cmd_arguments$fasta_file)
  gc()
  cat("================================================================================\n")
  if (log_messages != "") cat("03 / 14 \nFasta sizes: ", sapply(fasta_content, length), "\n", file = log_messages, append = TRUE)
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "03 Fasta loaded")
  times$data_type <- append(times$data_type, "Fasta total length")
  times$data_value <- append(times$data_value, sum(sapply(fasta_content, length)))

  ### 04 / 14 Calculate repeat scores for each sequence =================================================================
  cat(" 04 / 13 Calculating repeat scores for each sequence             ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  pb <- txtProgressBar(min = 0, max = length(fasta_content), style = 1)
  repeat_scores <- list()
  for (i in seq_along(fasta_content)) {
    repeat_scores <- append(repeat_scores, list(sequence_window_score(fasta_content[[i]], window_size, kmer))) # it's parallel
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
  }
  close(pb)
  gc()
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 04 sequence window score")
  times$data_type <- append(times$data_type, "Fasta total length")
  times$data_value <- append(times$data_value, sum(sapply(fasta_content, length)))

  ### 05 / 14 Identify regions with high repeat content and merge into a df =============================================
  cat(" 05 / 13 Identifying regions with high repeat content            ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  repetitive_regions <- data.frame(starts = NULL, ends = NULL, scores = NULL, seqID = NULL, numID = NULL)
  for (i in seq_along(repeat_scores)) {
    regions_of_sequence <- merge_windows(list_of_scores = repeat_scores[[i]], window_size = window_size, sequence_full_length = length(fasta_content[[i]]), log_messages)
    if (nrow(regions_of_sequence) != 0) {
      regions_of_sequence$seqID <- names(fasta_content)[[i]]
      regions_of_sequence$numID <- i
      repetitive_regions <- rbind(repetitive_regions, regions_of_sequence)
    }
  }
  write.csv(x = repetitive_regions, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_regarrays.csv")), row.names = FALSE)

  if (!inherits(repetitive_regions, "data.frame")) stop("No regions with repeats identified")
  if (nrow(repetitive_regions) == 0) stop("No regions with repeats identified")
  cat("================================================================================\n")
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 05 merge windows into regions")
  times$data_type <- append(times$data_type, "Number of windows to merge")
  times$data_value <- append(times$data_value, sum(sapply(repeat_scores, length)))

  ### 06 / 14 Split regions into arrays =================================================================================
  cat(" 06 / 13 Identifying individual arrays with repeats              ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  region_sizes <- repetitive_regions$ends - repetitive_regions$starts
  progress_values <- region_sizes / sum(region_sizes)
  pb <- txtProgressBar(min = 0, max = 1, style = 1)
  print("where error")
  arrays <- foreach (i = seq_len(nrow(repetitive_regions)),
                     .combine = rbind,
                     .export = c("split_and_check_arrays", "extract_kmers", "collapse_kmers", "genomic_bins_starts", "consensus_N", "write_align_read", "seq_win_score_int")) %dopar% {
    out <- split_and_check_arrays(start = repetitive_regions$starts[i],
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
                                  sink_output = FALSE,
                                  kmer = kmer)
    setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
    return(out)
  }
  print("where error")
  close(pb)
  gc()
  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_aregarrays.csv")), row.names = FALSE)


  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 06 split regions into arrays")
  times$data_type <- append(times$data_type, "Total length of arrays")
  times$data_value <- append(times$data_value, sum(arrays$end - arrays$start))

  ### 07 / 14 Shift representative repeats and apply templates ==========================================================
  cat(" 07 / 13 Shifting representative and comparing templates         ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  pb <- txtProgressBar(min = 0, max = nrow(arrays), style = 1)
  if (cmd_arguments$templates != 0) {
    templates <- read_fasta_and_list(cmd_arguments$templates)
    length_templates <- length(templates)
    if (length(templates) == 0) stop("No templates found within the template file")
  } else {
    templates <- 0
    length_templates <- 0
  }
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 07 get templates")
  times$data_type <- append(times$data_type, "Number of templates")
  times$data_value <- append(times$data_value, length_templates)

  arrays$representative <- foreach (i = seq_len(nrow(arrays)), .combine = c, .export = c("shift_and_compare", "shift_sequence", "compare_circular", "rev_comp_string", "kmer_hash_score")) %dopar% {
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    if (!inherits(arrays$representative[i], "character")) return("_")
    return(shift_and_compare(arrays$representative[i], templates))
  }
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 07 shift representatives and apply templates")
  times$data_type <- append(times$data_type, "Number of templates times nrow arrays")
  times$data_value <- append(times$data_value, length_templates * nrow(arrays))

  arrays$class <- ""
  for (i in seq_len(nrow(arrays))) {
    if (arrays$representative[i] == "_") {
      arrays$representative[i] <- ""
      next()
    }
    arrays$class[i] <- strsplit(arrays$representative[i], split = "_")[[1]][1]
    arrays$representative[i] <- strsplit(arrays$representative[i], split = "_")[[1]][2]
  }
  close(pb)
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 07")
  times$data_type <- append(times$data_type, "Arrays nrow")
  times$data_value <- append(times$data_value, nrow(arrays))

  ### 08 / 14 Classify unclassified and shift ===========================================================================
  cat(" 08 / 13 Classifying remaining representative repeats            ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  arrays <- classify_repeats(repeat_df = arrays)
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 08 classify repeats")
  times$data_type <- append(times$data_type, "Arrays nrow")
  times$data_value <- append(times$data_value, nrow(arrays))

  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_arrays.csv")), row.names = FALSE)

  classes <- unique(arrays$class)
  classes <- classes[!(classes %in% c(names(templates), "none_identified"))]
  if (length(classes) != 0) {
    pb <- txtProgressBar(style = 1, min = 0, max = length(classes))
    arrays_t <- foreach (i = seq_along(classes), .combine = rbind, .export = c("compare_kmer_grep", "shift_classes", "compare_circular", "rev_comp_string")) %dopar% {
      arrays_class <- arrays[arrays$class == classes[i], ]
      arrays_class$representative <- shift_classes(arrays_class, kmer = 6)
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      return(arrays_class)
    }
    arrays <- rbind(arrays_t, arrays[which(arrays$class %in% c(names(templates), "none_identified")), ])
    remove(arrays_t)
  } else {
    cat("================================================================================")
  }
  cat("\n")
  gc()

  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 08 shift classes")
  times$data_type <- append(times$data_type, "Unique classes number")
  times$data_value <- append(times$data_value, length(unique(arrays$class)))

  ### 09 / 14 Map repeats ===============================================================================================
  cat(" 09 / 13 Mapping array representatives                           ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  array_sizes <- arrays$end - arrays$start
  progress_values <- array_sizes / sum(array_sizes)
  pb <- txtProgressBar(style = 1)
  repeats <- foreach (i = seq_len(nrow(arrays)), .combine = rbind, .export = c("fill_gaps", "write_align_read", "consensus_N", "read_and_format_nhmmer", "handle_overlaps", "handle_gaps", "export_gff", "map_nhmmer", "map_default", "rev_comp_string")) %dopar% {
    time_report_df <- as.numeric(Sys.time())
    default_df = data.frame(seqID = vector(mode = "character"),
                        arrayID = vector(mode = "numeric"),
                        start = vector(mode = "numeric"),
                        end = vector(mode = "numeric"),
                        strand = vector(mode = "character"),
                        score = vector(mode = "numeric"),
                        eval = vector(mode = "numeric"),
                        width = vector(mode = "numeric"),
                        class = vector(mode = "character"),
                        representative = vector(mode = "character"),
                        score_template = vector(mode = "numeric"))
    if (arrays$representative[i] == "") {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      return(default_df)
    }
    if (arrays$top_N[i] >= 14) {
      ## nhmmer for repeats of 14+ bp =======================================
      repeats_df <- map_nhmmer(cmd_arguments$output_folder, i, arrays$representative[i], arrays$seqID[i], arrays$start[i],
                               arrays$end[i], fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], nhmmer_dir)
    } else {
      ## matchpattern for shorter ===========================================
      repeats_df <- map_default(i, arrays$representative[i], arrays$seqID[i], arrays$start[i], paste(fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], collapse = ""))
    }
    time_report_df <- c(time_report_df, as.numeric(Sys.time()))

    if (nrow(repeats_df) < 2) {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      return(default_df)
    }
    ## add width and class ============================================================
    repeats_df$width <- repeats_df$end - repeats_df$start + 1
    repeats_df$class <- arrays$class[i]
    ## Handle overlaps ======================================================

    repeats_df <- handle_overlaps(repeats_df, overlap_threshold = 0.1)
    if (nrow(repeats_df) < 2) {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      return(default_df)
    }
    time_report_df <- c(time_report_df, as.numeric(Sys.time()))
    ## handle gaps if proper array ==========================================
    repeats_df <- handle_gaps(repeats_df, representative_len = arrays$top_N[i])

    time_report_df <- c(time_report_df, as.numeric(Sys.time()))

    ## Return nothing if handle_gaps removed all repeats ====================
    if (nrow(repeats_df) < 2) {
      setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
      gc()
      return(default_df)
    }
    ## Recalculate representative ===========================================
    # TODO: Use more repeats for long arrays, this is stringent and good for most arrays, but some deserve a better recalculation
    max_repeats_to_align <- 15
    min_repeats_to_recalculate <- 10
    repeats_df$representative <- arrays$representative[i]
    repeats_df$score_template <- -1
    sample_IDs <- which(repeats_df$strand != ".")
    if (length(sample_IDs) >= min_repeats_to_recalculate) {
      if (length(sample_IDs) > max_repeats_to_align) sample_IDs <- sample(sample_IDs, max_repeats_to_align)
      repeats_seq <- unlist(lapply(sample_IDs, function(X) paste0(fasta_content[[arrays$numID[i]]][repeats_df$start[X] : repeats_df$end[X]], collapse = "")))
      strands <- repeats_df$strand[sample_IDs]
      repeats_seq[which(strands == "-")] <- unlist(lapply(repeats_seq[which(strands == "-")], rev_comp_string))
      alignment <- write_align_read(mafft_exe = mafft_dir,
                                    temp_dir = cmd_arguments$output_folder,
                                    sequences = repeats_seq,
                                    name = paste(basename(cmd_arguments$fasta_file), arrays$seqID[i], arrays$numID[i], i, runif(1, 0, 1), sep = "_"))
      consensus <- consensus_N(alignment, arrays$top_N[i])
      if (length(consensus) != 0) repeats_df$representative <- consensus
    }
    time_report_df <- c(time_report_df, as.numeric(Sys.time()))

    ## check if short gaps contain the repeat ===============================
    repeats_df <- fill_gaps(repeats_df, fasta_content[[arrays$numID[i]]][arrays$start[i] : arrays$end[i]], arrays$start[i])
    time_report_df <- c(time_report_df, as.numeric(Sys.time()))

    ## Change to edit distance based score ==================================
    # TODO: use only unique repeats to recalculate, should make it run faster for arrays with many repeats (where many are also identical)
    repeats_seq <- unlist(lapply(seq_len(nrow(repeats_df)), function(X) paste0(fasta_content[[arrays$numID[i]]][repeats_df$start[X] : repeats_df$end[X]], collapse = "")))
    costs <- list(insertions = 1, deletions = 1, substitutions = 1)
    if (sum(repeats_df$strand == "+") > 0) repeats_df$score[repeats_df$strand == "+"] <- adist(repeats_df$representative[1], repeats_seq[repeats_df$strand == "+"], costs)[1, ]  / nchar(repeats_df$representative[1]) * 100
    if (sum(repeats_df$strand == "-") > 0) repeats_df$score[repeats_df$strand == "-"] <- adist(rev_comp_string(repeats_df$representative[1]), repeats_seq[repeats_df$strand == "-"])[1, ]  / nchar(repeats_df$representative[1]) * 100
    if (arrays$class[i] %in% names(templates)) {
      template <- paste(templates[[which(arrays$class[i] == names(templates))]], collapse = "")
      if (sum(repeats_df$strand == "+") > 0) repeats_df$score_template[repeats_df$strand == "+"] <- adist(template, repeats_seq[repeats_df$strand == "+"], costs)[1, ]  / nchar(template) * 100
      if (sum(repeats_df$strand == "-") > 0) repeats_df$score_template[repeats_df$strand == "-"] <- adist(rev_comp_string(template), repeats_seq[repeats_df$strand == "-"])[1, ]  / nchar(template) * 100
    }
    ### Correct repeats split into two by nhmmer ============================
    score_min_to_merge <- 30
    size_max_to_merge <- 1.0
    i <- 1
    while(i < nrow(repeats_df)) {
      if (repeats_df$score[i] > score_min_to_merge &&
            repeats_df$score[i + 1] > score_min_to_merge &&
            (sum(repeats_df$width[i : (i + 1)]) < (size_max_to_merge * nchar(repeats_df$representative[1])))) {
        # both have high score and are short
        if ((repeats_df$strand[i] == "+") && (repeats_df$strand[i + 1] == "+")) {
          new_score = adist(repeats_df$representative[1], paste0(repeats_seq[i : (i + 1)], collapse = ""), costs)[1, ]  / nchar(repeats_df$representative[1]) * 100
          if (new_score < min(repeats_df$score[i : (i + 1)])) {
            repeats_df$end[i] = repeats_df$end[i + 1]
            repeats_df = repeats_df[-(i + 1),]
            repeats_seq = repeats_seq[-(i + 1)]
            repeats_seq[i] = paste0(fasta_content[[arrays$numID[i]]][repeats_df$start[i] : repeats_df$end[i]], collapse = "")
            repeats_df$score[i] = new_score
            if (arrays$class[i] %in% names(templates)) {
              repeats_df$score_template[i] <- adist(template, repeats_seq[i])[1, ]  / nchar(template) * 100
            }
          }
        } else if ((repeats_df$strand[i] == "-") && (repeats_df$strand[i + 1] == "-")) {
          new_score = adist(rev_comp_string(repeats_df$representative[1]), paste0(repeats_seq[i : (i + 1)], collapse = ""), costs)[1, ]  / nchar(repeats_df$representative[1]) * 100
          if (new_score < min(repeats_df$score[i : (i + 1)])) {
            repeats_df$end[i] = repeats_df$end[i + 1]
            repeats_df = repeats_df[-(i + 1), ]
            repeats_seq = repeats_seq[-(i + 1)]
            repeats_seq[i] = paste0(fasta_content[[arrays$numID[i]]][repeats_df$start[i] : repeats_df$end[i]], collapse = "")
            repeats_df$score[i] = new_score
            if (arrays$class[i] %in% names(templates)) {
              repeats_df$score_template[i] <- adist(rev_comp_string(template), repeats_seq[i])[1, ]  / nchar(template) * 100
            }
          }
        }
      }
      i = i + 1
    }
    time_report_df <- c(time_report_df, as.numeric(Sys.time()))
    gc()
    setTxtProgressBar(pb, getTxtProgressBar(pb) + progress_values[i])
    if (nrow(repeats_df) < 2) {
      return(default_df)
    }
    repeats_df <- repeats_df[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class", "representative", "score_template")]
    repeats_df$seqID <- as.character(repeats_df$seqID)
    repeats_df$arrayID <- as.numeric(repeats_df$arrayID)
    repeats_df$start <- as.numeric(repeats_df$start)
    repeats_df$end <- as.numeric(repeats_df$end)
    repeats_df$strand <- as.character(repeats_df$strand)
    repeats_df$score <- as.numeric(repeats_df$score)
    repeats_df$eval <- as.numeric(repeats_df$eval)
    repeats_df$width <- as.numeric(repeats_df$width)
    repeats_df$class <- as.character(repeats_df$class)
    repeats_df$representative <- as.character(repeats_df$representative)
    repeats_df$score_template <- as.numeric(repeats_df$score_template)
    return(repeats_df)
  }
  close(pb)
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 09 map repeats")
  times$data_type <- append(times$data_type, "Repeats number")
  times$data_value <- append(times$data_value, nrow(repeats))

  ### 10 / 14 ===========================================================================================================

  ### 11 / 14 Summarise array information ===============================================================================
  cat(" 11 / 13 Summarising array information                           ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")

  for (i in seq_len((nrow(arrays)))) {
    if (sum((repeats$arrayID == i) > 0)) {
      arrays$representative[i] <- repeats$representative[which(repeats$arrayID == i)[1]]
    }
  }
  repeats <- repeats[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class", "score_template")]

  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 11 reassign array representatives")
  times$data_type <- append(times$data_type, "Arrays nrow")
  times$data_value <- append(times$data_value, nrow(arrays))

  arrays$repeats_number <- 0
  arrays$median_repeat_width <- 0
  arrays$median_score <- -1

  for (i in seq_len(nrow(arrays))) {
    repeats_temp <- repeats[repeats$arrayID == i, ]
    if (nrow(repeats_temp) > 0) {
      arrays$repeats_number[i] <- nrow(repeats_temp)
      arrays$median_repeat_width[i] <- ceiling(median(repeats_temp$width))
      arrays$median_score[i] <- ceiling(median(repeats_temp$score))
    }
  }
  arrays <- arrays[arrays$repeats_number != 0, ]
  gc()
  cat("================================================================================\n")

  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 11 summarise array info")
  times$data_type <- append(times$data_type, "Arrays nrow")
  times$data_value <- append(times$data_value, nrow(arrays))

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
             attributes = c(9, 10, 11),
             attribute.names = c("Name=", "Repeat_no=", "Repeat_median_width="))
  cat("================================================================================\n")
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 12 saved array info")
  times$data_type <- append(times$data_type, "Arrays nrow")
  times$data_value <- append(times$data_value, nrow(arrays))

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
             attributes = c(9, 6, 10),
             attribute.names = c("Name=", "Arry_EDS=", "Family_EDS="))
  cat("================================================================================\n")
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 13 saved repeats info")
  times$data_type <- append(times$data_type, "Repeats number")
  times$data_value <- append(times$data_value, nrow(repeats))

  ### 14 / 14 Done ======================================================================================================
  if (report_runtime) {
    times$time_passed <- 0
    times$time_per_Mbp <- 0
    times$time_per_event_data_value <- 0
    for (i in 2 : length(times$time)) {
      times$time_passed <- append(times$time_passed, (times$time[i] - times$time[i - 1]))
      times$time_per_Mbp <- append(times$time_per_Mbp, (1000000 * times$time_passed[i] / times$data_value[2]))
      times$time_per_event_data_value <- append(times$time_per_event_data_value, (1000000 * times$time_passed[i] / times$data_value[i]))
    }
    write.csv(x = times, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_run_time.csv")), row.names = FALSE)
  }
  if (log_messages != "") cat("14 / 14 \n## Done ##\nTime:         ", date(), "\n", file = log_messages, append = TRUE)

  stopCluster(cl)
  gc()
}
