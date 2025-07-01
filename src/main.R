main <- function(cmd_arguments) {
  cat("################################################################################\n")
  cat("###           TRASH: workspace initialised                       ")
  cat(Sys.time())
  cat("  ###\n")
  cat("################################################################################\n")
  # TODO: remove this development settings
  if (Sys.info()["sysname"] == "Windows") {
    mafft_dir <- "../dep/mafft-7.520-win64-signed/mafft-win/mafft.bat"
    # nhmmer_dir <- "../dep/hmmer/nhmmer.exe"
    nhmmer_dir <- "C:/cygwin64/home/Piotr Włodzimierz/hmmer/hmmer-3.4/src/nhmmer.exe"
    # nhmmer_dir <- "C:/cygwin64/home/vlothec/bin/nhmmer.exe"
  } else {
    mafft_dir <- "mafft"
    nhmmer_dir <- "nhmmer"
  }
  cat("================================================================================\n")

  ### 01 / 14 Start workers =============================================================================================
  log_messages <- file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_main_log_file.txt")) #TODO make into a flag
  log_messages <- ""

  cl <- makeCluster(cmd_arguments$cores_no,
                    # outfile="",
                    homogeneous = TRUE)
  clusterEvalQ(cl, .libPaths(c(.libPaths(), gsub("src", "R_libs", getwd()))))
  clusterEvalQ(cl, sink())
  registerDoParallel(cl)
  # foreach::foreach (i = 1 : getDoParWorkers()) %dopar% {
  #   # set.seed(0) # Sets random seed for reproducibility
  #   # setwd(cmd_arguments$output_folder)
  #   # options(error = function() {traceback(2, max.lines=500); if(!interactive()) quit(save="no", status=1, runLast=T)})
  # }

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
  window_size <- round((cmd_arguments$max_rep_size + kmer) * 1.1)
  report_runtime <- TRUE
  add_sequence_info <- TRUE

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
  if(length(fasta_content) == 0) {warning("Fasta could no be read or is empty"); return(1)}

  ### 04 / 14 Calculate repeat scores for each sequence =================================================================
  cat(" 04 / 13 Calculating repeat scores for each sequence             ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  chromosome_lengths <- unlist(lapply(seq_along(fasta_content), function(X) length(fasta_content[[X]])))
  cat("  Assembly total length:\t", round(sum(chromosome_lengths) / 1000000, 1), "Mbp \n")
  cat("  Sequences count:\t\t\t", length(chromosome_lengths), " \n")
  cat("  Sequences names:\t\t\t", names(fasta_content), "\n")
  cat("  Sequences lengths (bp):\t", chromosome_lengths, "\n\n")

  if (length(names(fasta_content)) != length(unique(names(fasta_content)))) {
    warning(paste0("\nWARNING: Sequence names in the ", basename(cmd_arguments$fasta_file), " fasta file are not unique \n They were appended to avoid assignment errors \n"))
    cat(paste0("\nWARNING: Sequence names in the ", basename(cmd_arguments$fasta_file), " fasta file are not unique \n They were appended to avoid assignment errors \n"))
    cat("\n", "Adjustments made: \n", paste = "")

    for (i in seq_along(unique(names(fasta_content)))) {
      if (sum(names(fasta_content) %in% unique(names(fasta_content))[i]) > 1) {
        cat("Old names:",  names(fasta_content)[names(fasta_content) == unique(names(fasta_content))[i]])
        cat("\n New names:",   paste0(unique(names(fasta_content))[i], 1:sum(names(fasta_content) %in% unique(names(fasta_content))[i])), "\n\n\n")
        names(fasta_content)[names(fasta_content) == unique(names(fasta_content))[i]] <- paste0(unique(names(fasta_content))[i], 1:sum(names(fasta_content) %in% unique(names(fasta_content))[i]))
      }
    }
  }

  repeat_scores <- list()
  for (i in seq_along(fasta_content)) {
    cat("  Fasta sequence ", i, ": ", names(fasta_content)[i], " \t", sep = "")
    repeat_scores <- append(repeat_scores, list(sequence_window_score(fasta_content[[i]], window_size, kmer, output_dir = cmd_arguments$output_folder))) # it's parallel inside
  }
  cat("================================================================================\n")
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
    if (length(repeat_scores[[i]]) == 0) next
    regions_of_sequence <- merge_windows(list_of_scores = repeat_scores[[i]], window_size = window_size, sequence_full_length = length(fasta_content[[i]]), log_messages)
    if (nrow(regions_of_sequence) != 0) {
      regions_of_sequence$seqID <- names(fasta_content)[[i]]
      regions_of_sequence$numID <- i
      repetitive_regions <- rbind(repetitive_regions, regions_of_sequence)
    }
  }
  write.csv(x = repetitive_regions, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_regarrays.csv")), row.names = FALSE)

  if (!inherits(repetitive_regions, "data.frame")) {
    print("No regions with repeats identified")
    return(0)
  }
  if (nrow(repetitive_regions) == 0) {
    print("No regions with repeats identified")
    return(0)
  }
  cat("================================================================================\n")
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 05 merge windows into regions")
  times$data_type <- append(times$data_type, "Number of windows to merge")
  times$data_value <- append(times$data_value, sum(sapply(repeat_scores, length)))
  remove(repeat_scores)
  gc()
  # warnings()
  ### 06 / 14 Split regions into arrays =================================================================================
  cat(" 06 / 13 Identifying individual arrays with repeats              ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  date <- Sys.Date()
  regions_per_chunk <- 100
  # region_sizes <- repetitive_regions$ends - repetitive_regions$starts
  arrays <- NULL
  for (i in seq_along(fasta_content)) {
    cat("  Fasta sequence ", i, ": ", names(fasta_content)[i], " \t", sep = "")
    repetitive_regions_chr <- repetitive_regions[repetitive_regions$numID == i,]
    cat("Repetitive regions in the sequence: ", length(repetitive_regions_chr), " \t", sep = "")
    if(nrow(repetitive_regions_chr) == 0) {
      cat("\n")
      next
    }
    region_chunk <- seq(1, nrow(repetitive_regions_chr), regions_per_chunk) #divide into up to 100 data frame entries chunks on each chromosome, so up to 100 parallel, too much of a fasta is not good sent into the parallel
    cat("Chunks to complete: ", length(region_chunk), ". Finished: ", sep = "")
    region_chunk <- c(region_chunk, (nrow(repetitive_regions_chr) + 1))
    for(j in 1 : (length(region_chunk) - 1)) {
      sequence_substring <- fasta_content[[i]][repetitive_regions_chr$starts[region_chunk[j]] : (repetitive_regions_chr$ends[region_chunk[j + 1] - 1])]
      start_adjust <- repetitive_regions_chr$starts[region_chunk[j]] - 1
      foreach::foreach (k = (region_chunk[j] : (region_chunk[j+1] - 1)),
                    .export = c("split_and_check_arrays", "extract_kmers", "collapse_kmers", "genomic_bins_starts", "consensus_N", "write_align_read", "seq_win_score_int")) %dopar% {
        
        out <- split_and_check_arrays(start = repetitive_regions_chr$starts[k],
                                      end = repetitive_regions_chr$ends[k],
                                      sequence = sequence_substring[(repetitive_regions_chr$starts[k] - start_adjust) : (repetitive_regions_chr$ends[k] - start_adjust)],
                                      seqID = repetitive_regions_chr$seqID[k],
                                      numID = repetitive_regions_chr$numID[k],
                                      arrID = k,
                                      max_repeat = cmd_arguments$max_rep_size,
                                      min_repeat = cmd_arguments$min_rep_size,
                                      mafft = mafft_dir,
                                      temp_dir = cmd_arguments$output_folder,
                                      src_dir = getwd(),
                                      sink_output = FALSE,
                                      kmer = kmer)
        save(out, file = paste0(cmd_arguments$output_folder, "/", i, "_", j, "_", k, "_", date, "_06_data"))
        remove(out)
        gc()
        # warnings()
        # return(out)
      }
      cat(j, "")
      for (k in (region_chunk[j] : (region_chunk[j+1] - 1))) {
        load(paste0(cmd_arguments$output_folder, "/", i, "_", j, "_", k, "_", date, "_06_data"))
        unlink(paste0(cmd_arguments$output_folder, "/", i, "_", j, "_", k, "_", date, "_06_data"))
        arrays <- rbind(arrays, out)
        remove(out)
      }
    }
    cat("\n")
  }
  remove(repetitive_regions)
  gc()
  cat("================================================================================\n")
  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_aregarrays.csv")), row.names = FALSE)


  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 06 split regions into arrays")
  times$data_type <- append(times$data_type, "Total length of arrays")
  times$data_value <- append(times$data_value, sum(arrays$end - arrays$start))
  # warnings()
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

  arrays$representative <- foreach::foreach (i = seq_len(nrow(arrays)),
                                    .combine = c,
                                    .export = c("shift_and_compare", "shift_sequence", "compare_circular", "rev_comp_string", "kmer_hash_score")) %dopar% {
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
    if (arrays$representative[i] == "_split_") {
      arrays$representative[i] <- ""
      next()
    }
    arrays$class[i] <- strsplit(arrays$representative[i], split = "_split_")[[1]][1]
    arrays$representative[i] <- strsplit(arrays$representative[i], split = "_split_")[[1]][2]
  }
  close(pb)
  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 07")
  times$data_type <- append(times$data_type, "Arrays nrow")
  times$data_value <- append(times$data_value, nrow(arrays))
  # warnings()
  ### 08 / 14 Classify unclassified and shift ===========================================================================
  cat(" 08 / 13 Classifying remaining representative repeats            ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  date <- Sys.Date()
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
    # arrays_t <- foreach::foreach (i = seq_along(classes), .combine = rbind, .export = c("compare_kmer_grep", "shift_classes", "compare_circular", "rev_comp_string")) %dopar% {
    foreach::foreach (i = seq_along(classes), .export = c("compare_kmer_grep", "shift_classes", "compare_circular", "rev_comp_string")) %dopar% {
      arrays_class <- arrays[arrays$class == classes[i], ]
      arrays_class$representative <- shift_classes(arrays_class, kmer = 6)
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      save(arrays_class, file = paste0(cmd_arguments$output_folder, "/", i, "_", date, "_08_data"))
      remove(arrays_class)
      gc()
      # return(arrays_class)
    }
    arrays_t <- NULL
    for (i in seq_along(classes)) {
      load(paste0(cmd_arguments$output_folder, "/", i, "_", date, "_08_data"))
      unlink(paste0(cmd_arguments$output_folder, "/", i, "_", date, "_08_data"))
      arrays_t <- rbind(arrays_t, arrays_class)
      remove(arrays_class)
    }
    arrays <- rbind(arrays_t, arrays[which(arrays$class %in% c(names(templates), "none_identified")), ])
    close(pb)
    remove(arrays_t)
    gc()
  } else {
    cat("================================================================================\n")
  }
  arrays <- arrays[order(arrays$start), ]
  arrays <- arrays[order(arrays$seqID), ]
  arrays$array_num_ID <- seq_len(nrow(arrays))

  arrays_no_representative <- arrays[arrays$class == "none_identified", ]
  arrays <- arrays[arrays$class != "none_identified", ]
  write.csv(x = arrays_no_representative, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_no_repeats_arrays.csv")), row.names = FALSE)
  write.csv(x = arrays, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_classarrays.csv")), row.names = FALSE)

  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 08 shift classes")
  times$data_type <- append(times$data_type, "Unique classes number")
  times$data_value <- append(times$data_value, length(unique(arrays$class)))

  remove(arrays_no_representative)
  gc()

  if(nrow(arrays) == 0) {
    cat("No arrays with tandem repeats found under the settings\n")
    return(0)
  }

  ### 09 / 14 Map repeats ===============================================================================================
  cat(" 09 / 13 Mapping array representatives                           ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  
  repeats <- NULL
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
  arrays_per_chunk <- 100

  for(chromosome in seq_along(fasta_content)) {
    ## For each chromosome ========================================================
    cat("  Fasta sequence ", chromosome, ": ", names(fasta_content)[chromosome], " \t", sep = "")
    arrays_chr <- arrays[arrays$numID == chromosome,]
    cat("Arrays in the sequence: ", nrow(arrays_chr), " \t", sep = "")
    if(nrow(arrays_chr) == 0) {
      cat("\n")
      next
    }
    array_chunk <- seq(1, nrow(arrays_chr), arrays_per_chunk) #divide into up to 100 data frame entries chunks on each chromosome, so up to 100 parallel, too much of a fasta is not good sent into the parallel
    cat("Chunks to complete: ", length(array_chunk), ". Finished: ", sep = "")
    array_chunk <- c(array_chunk, (nrow(arrays_chr) + 1))
    for(j in 1 : (length(array_chunk) - 1)) {
      sequence_substring <- fasta_content[[chromosome]][arrays_chr$start[array_chunk[j]] : (arrays_chr$end[array_chunk[j + 1] - 1])]
      adjust_start <- arrays_chr$start[array_chunk[j]] - 1
      arrays_chunk_IDs <- array_chunk[j] : (array_chunk[j + 1] - 1)
      foreach::foreach (i = arrays_chunk_IDs,
                        .combine = rbind,
                        # .packages = c("Biostrings", "seqinr", "msa"),
                        .export = c("find_edge_best_start_end", "handle_edge_repeat", "fill_gaps", "write_align_read", "consensus_N", "read_and_format_nhmmer", "handle_overlaps", "handle_gaps", "export_gff", "map_nhmmer", "map_default", "rev_comp_string")) %dopar% {
        if (arrays_chr$representative[i] == "") {
          cat(i, "")
          return(0)
        }
        array_sequence <- sequence_substring[(arrays_chr$start[i] - adjust_start) : (arrays_chr$end[i] - adjust_start)]
        cat(i, "_ ", sep = "")
        if (arrays_chr$top_N[i] >= 14) {
          # nhmmer for repeats of 14+ bp =======================================
          repeats_df <- map_nhmmer(cmd_arguments$output_folder, arrayID = arrays_chr$array_num_ID[i], arrays_chr$representative[i], arrays_chr$seqID[i], arrays_chr$start[i],
                                   arrays_chr$end[i], array_sequence, nhmmer_dir)
        } else {
          # matchpattern for shorter ===========================================
          repeats_df <- map_default(arrayID = arrays_chr$array_num_ID[i], arrays_chr$representative[i], arrays_chr$seqID[i], arrays_chr$start[i], paste(array_sequence, collapse = ""))
        }
        if (nrow(repeats_df) < 2) {
          cat(i, "")
          return(0)
        }
        # add width and class ============================================================
        repeats_df$width <- repeats_df$end - repeats_df$start + 1
        repeats_df$class <- arrays_chr$class[i]
        # Handle overlaps ======================================================
        repeats_df <- handle_overlaps(repeats_df, overlap_threshold = 0.1)
        if (nrow(repeats_df) < 3) {
          cat(i, "")
          return(0)
        }
        # handle gaps if proper array ==========================================
        repeats_df <- handle_gaps(repeats_df, representative_len = arrays_chr$top_N[i])
        # Return nothing if handle_gaps removed all repeats ====================
        if (nrow(repeats_df) < 3) {
          cat(i, "")
          return(0)
        }
        # Recalculate representative ===========================================
        # TODO: Use more repeats for long arrays, this is stringent and good for most arrays, but some deserve a better recalculation
        max_repeats_to_align <- 15
        min_repeats_to_recalculate <- 10
        repeats_df$representative <- arrays_chr$representative[i]
        repeats_df$score_template <- -1
        sample_IDs <- which(repeats_df$strand != ".")
        if (length(sample_IDs) >= min_repeats_to_recalculate) {
          if (length(sample_IDs) > max_repeats_to_align) sample_IDs <- sample(sample_IDs, max_repeats_to_align)
          repeats_seq <- unlist(lapply(sample_IDs, function(X) paste0(sequence_substring[(repeats_df$start[X] - adjust_start) : (repeats_df$end[X] - adjust_start)], collapse = "")))
          strands <- repeats_df$strand[sample_IDs]
          repeats_seq[which(strands == "-")] <- unlist(lapply(repeats_seq[which(strands == "-")], rev_comp_string))
          alignment <- write_align_read(mafft_exe = mafft_dir,
                                        temp_dir = cmd_arguments$output_folder,
                                        sequences = repeats_seq,
                                        name = paste(basename(cmd_arguments$fasta_file), arrays_chr$seqID[i], arrays_chr$array_num_ID[i], i, runif(1, 0, 1), sep = "_"))
          consensus <- consensus_N(alignment, arrays_chr$top_N[i])
          if (length(consensus) != 0) repeats_df$representative <- consensus
          remove(alignment, consensus, repeats_seq, strands)
          gc()
        }
        # check if short gaps contain the repeat ===============================
        repeats_df <- fill_gaps(repeats_df, array_sequence, arrays_chr$start[i])
        # Change to edit distance based score ==================================
        # TODO: use only unique repeats to recalculate, should make it run faster for arrays with many repeats (where many are also identical)
        repeats_seq <- unlist(lapply(seq_len(nrow(repeats_df)), function(X) paste0(sequence_substring[(repeats_df$start[X] - adjust_start) : (repeats_df$end[X] - adjust_start)], collapse = "")))
        costs <- list(insertions = 1, deletions = 1, substitutions = 1)
        if (sum(repeats_df$strand == "+") > 0) repeats_df$score[repeats_df$strand == "+"] <- adist(repeats_df$representative[1], repeats_seq[repeats_df$strand == "+"], costs)[1, ]  / nchar(repeats_df$representative[1]) * 100
        if (sum(repeats_df$strand == "-") > 0) repeats_df$score[repeats_df$strand == "-"] <- adist(rev_comp_string(repeats_df$representative[1]), repeats_seq[repeats_df$strand == "-"])[1, ]  / nchar(repeats_df$representative[1]) * 100
        if (arrays_chr$class[i] %in% names(templates)) {
          template <- paste(templates[[which(names(templates) == arrays_chr$class[i])]], collapse = "")
          if (sum(repeats_df$strand == "+") > 0) repeats_df$score_template[repeats_df$strand == "+"] <- adist(template, repeats_seq[repeats_df$strand == "+"], costs)[1, ]  / nchar(template) * 100
          if (sum(repeats_df$strand == "-") > 0) repeats_df$score_template[repeats_df$strand == "-"] <- adist(rev_comp_string(template), repeats_seq[repeats_df$strand == "-"])[1, ]  / nchar(template) * 100
        }
        # Correct repeats split into two by nhmmer ============================
        score_min_to_merge <- 30
        size_max_to_merge <- 1.0
        i_r <- 1
        while(i_r < nrow(repeats_df)) {
          if (repeats_df$score[i_r] > score_min_to_merge &&
              repeats_df$score[i_r + 1] > score_min_to_merge &&
              (sum(repeats_df$width[i_r : (i_r + 1)]) < (size_max_to_merge * nchar(repeats_df$representative[1])))) {
            # both have high score and are short
            if ((repeats_df$strand[i_r] == "+") && (repeats_df$strand[i_r + 1] == "+")) {
              new_score = adist(repeats_df$representative[1], paste0(repeats_seq[i_r : (i_r + 1)], collapse = ""), costs)[1, ]  / nchar(repeats_df$representative[1]) * 100
              if (new_score < min(repeats_df$score[i_r : (i_r + 1)])) {
                repeats_df$end[i_r] = repeats_df$end[i_r + 1]
                repeats_df = repeats_df[-(i_r + 1),]
                repeats_seq = repeats_seq[-(i_r + 1)]
                repeats_seq[i_r] = paste0(sequence_substring[(repeats_df$start[i_r] - adjust_start) : (repeats_df$end[i_r] - adjust_start)], collapse = "")
                repeats_df$score[i_r] = new_score
                if (arrays_chr$class[i_r] %in% names(templates)) {
                  template <- paste(templates[[which(names(templates) == arrays_chr$class[i_r])]], collapse = "")
                  repeats_df$score_template[i_r] <- adist(template, repeats_seq[i_r])[1, ]  / nchar(template) * 100
                }
              }
            } else if ((repeats_df$strand[i_r] == "-") && (repeats_df$strand[i_r + 1] == "-")) {
              new_score = adist(rev_comp_string(repeats_df$representative[1]), paste0(repeats_seq[i_r : (i_r + 1)], collapse = ""), costs)[1, ]  / nchar(repeats_df$representative[1]) * 100
              if (new_score < min(repeats_df$score[i_r : (i_r + 1)])) {
                repeats_df$end[i_r] = repeats_df$end[i_r + 1]
                repeats_df = repeats_df[-(i_r + 1), ]
                repeats_seq = repeats_seq[-(i_r + 1)]
                repeats_seq[i_r] = paste0(sequence_substring[(repeats_df$start[i_r] - adjust_start) : (repeats_df$end[i_r] - adjust_start)], collapse = "")
                repeats_df$score[i_r] = new_score
                if (arrays_chr$class[i_r] %in% names(templates)) {
                  template <- paste(templates[[which(names(templates) == arrays_chr$class[i_r])]], collapse = "")
                  repeats_df$score_template[i_r] <- adist(rev_comp_string(template), repeats_seq[i_r])[1, ]  / nchar(template) * 100
                }
              }
            }
          }
          i_r = i_r + 1
        }
        # Handle edge repeats =================================================
        repeats_df <- handle_edge_repeat(repeats_df, sequence_substring, adjust_start)
        # TODO: make sure the edge repeats have their template score recalculated too

        remove(repeats_seq)
        gc()
        if (nrow(repeats_df) < 3) {
          cat(i, "")
          return(0)
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
        save(repeats_df, file = paste0(cmd_arguments$output_folder, "/", i, "_", date, "_09_data"))
        remove(repeats_df)
        cat(i, "")
        return(0)
        gc()
      }
      for(i in arrays_chunk_IDs) {
        if(file.exists(paste0(cmd_arguments$output_folder, "/", i, "_", date, "_09_data"))) {
          load(paste0(cmd_arguments$output_folder, "/", i, "_", date, "_09_data"))
          unlink(paste0(cmd_arguments$output_folder, "/", i, "_", date, "_09_data"))
          if(sum(names(repeats) != names(repeats_df)) > 0) {
          print(paste(i, names(repeats), names(repeats_df)))
          print(str(repeats))
          print(str(repeats_df))
        }
        if(nrow(repeats_df) > 0) repeats <- rbind(repeats, repeats_df)
        remove(repeats_df)
        }
      }
    }
    cat("\n")
  }

  times$time <- append(times$time, as.numeric(Sys.time()))
  times$event <- append(times$event, "Finished 09 map repeats")
  times$data_type <- append(times$data_type, "Repeats number")
  times$data_value <- append(times$data_value, nrow(repeats))
  # warnings()
  ### 10 / 14 ===========================================================================================================

  ### 11 / 14 Summarise array information ===============================================================================
  cat(" 11 / 13 Summarising array information                           ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")

  for (i in seq_len(nrow(arrays))) {
    if (sum((repeats$arrayID == arrays$array_num_ID[i]) > 0)) {
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
    repeats_temp <- repeats[repeats$arrayID == arrays$array_num_ID[i], ]
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

  if (add_sequence_info) {
    repeats$sequence <- ""

    repeats$sequence <- unlist(lapply(seq_len(nrow(repeats)), function(X) paste0(fasta_content[[which(names(fasta_content) == repeats$seqID[X])]][repeats$start[X] : repeats$end[X]], collapse = "")))

    repeats$sequence[which(repeats$strand == "-")] <- unlist(lapply(repeats$sequence[which(repeats$strand == "-")], rev_comp_string))
    write.csv(x = repeats, file = file.path(cmd_arguments$output_folder, paste0(basename(cmd_arguments$fasta_file), "_repeats_with_seq.csv")), row.names = FALSE)
  }

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
