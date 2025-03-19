
# TODO: this is barrely changed hort script, adjusted for running analysis between chromosomes, many things (especially at the end) might not work as intended


horb <- function(cmd_arguments) {
  cat("TRASH HOR module: workspace initialised                        ")
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  # TODO: remove this development settings
  if (Sys.info()["sysname"] == "Windows") {
    mafft_dir <- normalizePath(file.path(gsub("src", "", getwd()), "dep", "mafft-7.520-win64-signed", "mafft-win", "mafft.bat"))
    hor_path <- normalizePath(file.path(gsub("src", "", getwd()), "dep", "new.hor"))
  } else {
    mafft_dir <- "mafft"
    hor_path <- normalizePath(file.path(gsub("src", "", getwd()), "dep", "HOR.V3.3"))
  }
  print(cmd_arguments)
  cat("================================================================================\n")

  ### 00 / 00 Start workers =============================================================================================

  set.seed(0) # Sets random seed for reproducibility

  ### 00 / 00 Settings ==================================================================================================
  report_runtime <- TRUE
  filter_different_blocks <- FALSE
  max_block_dif <- 1000
  threshold_SNV <- cmd_arguments$hor_threshold
  method_hor <- 2
  hor_c_verbose = FALSE
  SNV_per_kbp_in_red = 25 # from which hor$SNV_per_kbp value plotted lines will be red

  ### 00 / 00 Read repeats ==============================================================================================
  cat(paste0(" 00 / 00 ", paste(rep(" ", 50), collapse = "")))
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")

  repeats <- read.csv(cmd_arguments$repeats)
  repeatsB <- read.csv(cmd_arguments$repeatsB)

  cat("================================================================================\n")

  ### 03 / 00 Filter repeats and get sequences ==========================================================================
  cat(paste0(" 03 / 00 ", paste(rep(" ", 50), collapse = "")))
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")

  # seqID	arrayID	start	end	strand	score	eval	width	class	score_template          # repeats df col names


  repeats <- repeats[repeats$class == cmd_arguments$class, ]
  if (nrow(repeats) == 0) stop(paste0("No repeats of the class ", cmd_arguments$class, " found"))
  repeats <- repeats[repeats$seqID == cmd_arguments$chrA, ]
  if (nrow(repeats) == 0) stop(paste0("No repeats on the sequence ", cmd_arguments$chrA, " found"))
  repeats$strand[repeats$strand == "+"] <- "1"
  repeats$strand[repeats$strand == "-"] <- "2"


  repeatsB <- repeatsB[repeatsB$class == cmd_arguments$classB, ]
  if (nrow(repeatsB) == 0) stop(paste0("No repeatsB of the class ", cmd_arguments$classB, " found"))
  repeatsB <- repeatsB[repeatsB$seqID == cmd_arguments$chrB, ]
  if (nrow(repeatsB) == 0) stop(paste0("No repeatsB on the sequence ", cmd_arguments$chrB, " found"))
  repeatsB$strand[repeatsB$strand == "+"] <- "1"
  repeatsB$strand[repeatsB$strand == "-"] <- "2"

  split_after <- nrow(repeats)

  repeats <- rbind(repeats, repeatsB)
  
  threshold_SNV <- floor(threshold_SNV * median(repeats$width) / 100)
  cat("================================================================================\n")

  ### 04 / 00 Align all repeats =========================================================================================
  cat(paste0(" 04 / 00 ", paste(rep(" ", 50), collapse = "")))
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")

  alignment_name <- paste(cmd_arguments$class, cmd_arguments$classB, cmd_arguments$chrA, cmd_arguments$chrB, cmd_arguments$genomeA, cmd_arguments$genomeB, sep = "_")

  write_align_read(mafft_exe = mafft_dir,
                   temp_dir = cmd_arguments$output_folder,
                   sequences = repeats$sequence,
                   name = alignment_name,
                   seq_names = as.list(paste(seq_along(repeats$strand), "_D", repeats$strand, sep = "")),
                   use_mafft = TRUE,
                   remove_ali_mafft = FALSE,
                   return_ali_mafft = FALSE)

  alignment_file <- normalizePath(file.path(cmd_arguments$output_folder, paste0(alignment_name, "temp.aligned.fasta")))

  if (!file.exists(alignment_file)) stop("cannot find alignment file for HORs, stopping HOR calculation")
  if (file.size(alignment_file) == 0) stop("output alignment file for HORs empty, stopping HOR calculation")
  cat("================================================================================\n")

  ### 05 / 00 Run HORs ==================================================================================================
  cat(paste0(" 05 / 00 ", paste(rep(" ", 50), collapse = "")))
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")

  file_sep <- "/"
  sys_type <- "sh"
  if (.Platform$OS.type == "windows") {
    file_sep <- "\\"
    sys_type <- "cmd"
  }

  dir_0 <- paste0(cmd_arguments$output_folder, file_sep, collapse = "")

  run_arguments <- paste(shQuote(dir_0, type = sys_type),
                         shQuote(paste0(alignment_name, "temp.aligned.fasta"), type = sys_type),
                         split_after,
                         threshold_SNV,
                         cmd_arguments$hor_min_len,
                         method_hor, # 1 is all analysed, 2 is chrA vs chrB using split_after
                         sep = " ")

  print("Aligned, running HORs with arguments: ")
  print(hor_path)
  print(run_arguments)
  stdoutc <- NULL
  if (hor_c_verbose) stdoutc <- TRUE
  system2(command = hor_path,
          args = run_arguments, stdout = stdoutc, stderr = "")

  hor.output.file <- paste(cmd_arguments$output_folder, "/", "HORs_method_", method_hor, "_", paste0(alignment_name, "temp.aligned.fasta"), "_t_", threshold_SNV, "_c_", cmd_arguments$hor_min_len, ".csv", sep = "")

  hors <- read.csv(file = hor.output.file, header = TRUE)
  print(paste0("HORs identified: ", nrow(hors)))
  if (nrow(hors) > 1) {
    hors <- as.data.frame(hors)

    hors$start.A.bp <- repeats$start[hors$start_A]
    hors$start.B.bp <- repeats$start[hors$start_B]
    hors$end.A.bp <- repeats$end[hors$end_A]
    hors$end.B.bp <- repeats$end[hors$end_B]
    hors$chrA <- repeats$seqID[hors$start_A]
    hors$chrB <- repeats$seqID[hors$start_B]
    hors$block.size.in.units <- hors$end_A - hors$start_A + 1
    hors$block.A.size.bp <- hors$end.A.bp - hors$start.A.bp + 1
    hors$block.B.size.bp <- hors$end.B.bp - hors$start.B.bp + 1
    hors$SNV_per_kbp <- 1000 * hors$total_variant / ((hors$block.A.size.bp + hors$block.A.size.bp) / 2)
    if (filter_different_blocks) {
      for (p in rev(seq_len(nrow(hors)))) {
        block1size <- abs(hors$end.A.bp[p] - hors$start.A.bp[p]) + 1
        block2size <- abs(hors$end.B.bp[p] - hors$start.B.bp[p]) + 1
        if (abs(block2size - block1size) > max_block_dif) {
          hors <- hors[-p, ]
          p <- p - 1
        }
      }
    }
    write.csv(x = hors, file = file.path(cmd_arguments$output_folder, paste0("HORs_", cmd_arguments$class, "_", cmd_arguments$classB, "_", cmd_arguments$chrA, "_", cmd_arguments$chrB, "_", cmd_arguments$genomeA, "_", cmd_arguments$genomeB, ".csv", collapse = "")))

  } else {
    print("No HORs identified")
    file.remove(hor.output.file)
    file.remove(file.path(cmd_arguments$output_folder, paste0(alignment_name, "temp.aligned.fasta")))
    file.remove(paste(cmd_arguments$output_folder, "/", "log_", paste0(alignment_name, "temp.aligned.fasta"), "_t_", threshold_SNV, "_c_", cmd_arguments$hor_min_len, ".txt", sep = ""))
    summary <- data.frame(genomeA = cmd_arguments$genomeA,
                          genomeB = cmd_arguments$genomeB,
                          chrA = cmd_arguments$chrA,
                          chrB = cmd_arguments$chrB,
                          repA = cmd_arguments$class,
                          repB = cmd_arguments$classB,
                          repA_total = split_after,
                          repB_total = nrow(repeats) - split_after,
                          A_repeats_found_in_B = 0,
                          B_repeats_found_in_A = 0)

   write.csv(x = summary, file = file.path(cmd_arguments$output_folder, paste0("summary_of_hors_", cmd_arguments$class, "_", cmd_arguments$classB, "_", cmd_arguments$chrA, "_", cmd_arguments$chrB, "_", cmd_arguments$genomeA, "_", cmd_arguments$genomeB, ".csv", collapse = "")), row.names = FALSE)
  
    return(0)
  }
  file.remove(hor.output.file)
  file.remove(file.path(cmd_arguments$output_folder, paste0(alignment_name, "temp.aligned.fasta")))
  file.remove(paste(cmd_arguments$output_folder, "/", "log_", paste0(alignment_name, "temp.aligned.fasta"), "_t_", threshold_SNV, "_c_", cmd_arguments$hor_min_len, ".txt", sep = ""))
  cat("================================================================================\n")

  ### 06 / 00 Calculate bins and adjust starts ==========================================================================
  cat(paste0(" 06 / 00 ", paste(rep(" ", 50), collapse = "")))
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")

  if (cmd_arguments$plot_simple == "y") {
    unit.name <- "repeat ID"

    png(filename = file.path(cmd_arguments$output_folder, paste0("HORs_lines_simple_", cmd_arguments$class, "_", cmd_arguments$classB, "_", cmd_arguments$chrA, "_", cmd_arguments$chrB, "_", cmd_arguments$genomeA, "_", cmd_arguments$genomeB, ".png", collapse = "")), width = 6000, height = 6000, pointsize = 80)

    par(mar = c(4, 4, 4, 4), oma = c(1, 1, 1, 1))

    ax.len <- max(c(split_after, nrow(repeats) - split_after))

    hors = hors[order(hors$SNV_per_kbp, decreasing = TRUE), ]
    plot(x = NULL, y = NULL,
         xlab = paste0(cmd_arguments$chrA, ", ", unit.name),
         ylab = paste0(cmd_arguments$chrA, ", ", unit.name),
         xlim = c(0, ax.len),
         ylim = c(0, ax.len),
         pch = 19, cex = 0.1,
         main = paste0("HORs", cmd_arguments$class, " ", cmd_arguments$chrA))
    colours_SNV <- colorRampPalette(c("green", "yellow", "red"))(length(hors$SNV_per_kbp)) [findInterval(hors$SNV_per_kbp, seq(0, SNV_per_kbp_in_red, length.out = length(hors$SNV_per_kbp)))]
    while (TRUE) {
      if (ax.len >= 100000) {lwd_plot <- 1; break}
      if (ax.len >= 50000) {lwd_plot <- 2; break}
      if (ax.len >= 25000) {lwd_plot <- 3; break}
      if (ax.len >= 12500) {lwd_plot <- 4; break}
      lwd_plot <- 5
      break
    }
    for (j in seq_len(nrow(hors))) {
      lines(x = c(hors$start_A[j], hors$end_A[j]),
            y = c(hors$start_B[j] - split_after, hors$end_B[j] - split_after),
            pch = 19, lwd = lwd_plot, col = colours_SNV[j])
    }
    abline(h = nrow(repeats) - split_after, v = split_after)
    dev.off()



  } else {
    windows.to.check.if.any.repetas.are.in.and.plot <- 100000 # 100 Kbp
    ticks.every <- 100000
    unit.val <- 1000000
    unit.name <- "Mbp"

    #calculate plot regions
    bins.df <- calc_plot_regs(repeats$start, windows.to.check.if.any.repetas.are.in.and.plot)

    # Adjust repeats and HORs starts

    repeats$start.adjusted = repeats$start
    for (i in seq_len(nrow(bins.df))) {
      repeats$start.adjusted[repeats$start >= bins.df$bins.starts[i] & repeats$start < bins.df$bins.ends[i]] =
        repeats$start.adjusted[repeats$start >= bins.df$bins.starts[i] & repeats$start < bins.df$bins.ends[i]] - bins.df$correction[i]
      hors$start.A.bp[hors$start.A.bp >= bins.df$bins.starts[i] & hors$start.A.bp < bins.df$bins.ends[i]] =
        hors$start.A.bp[hors$start.A.bp >= bins.df$bins.starts[i] & hors$start.A.bp < bins.df$bins.ends[i]] - bins.df$correction[i]
      hors$end.A.bp[hors$end.A.bp >= bins.df$bins.starts[i] & hors$end.A.bp < bins.df$bins.ends[i]] =
        hors$end.A.bp[hors$end.A.bp >= bins.df$bins.starts[i] & hors$end.A.bp < bins.df$bins.ends[i]] - bins.df$correction[i]
      hors$start.B.bp[hors$start.B.bp >= bins.df$bins.starts[i] & hors$start.B.bp < bins.df$bins.ends[i]] =
        hors$start.B.bp[hors$start.B.bp >= bins.df$bins.starts[i] & hors$start.B.bp < bins.df$bins.ends[i]] - bins.df$correction[i]
      hors$end.B.bp[hors$end.B.bp >= bins.df$bins.starts[i] & hors$end.B.bp < bins.df$bins.ends[i]] =
        hors$end.B.bp[hors$end.B.bp >= bins.df$bins.starts[i] & hors$end.B.bp < bins.df$bins.ends[i]] - bins.df$correction[i]
    }

    ### 07 / 00 Make line plot ============================================================================================
    cat(paste0(" 07 / 00 ", paste(rep(" ", 50), collapse = "")))
    cat(Sys.time())
    cat("\n")
    cat("################################################################################\n")

    x.ax.len <- sum(bins.df$bin.size)
    y.ax.len <- x.ax.len
    ax.len <- max(c(x.ax.len, y.ax.len))

    png(filename = file.path(cmd_arguments$output_folder, paste0("HORs_lines_", cmd_arguments$class, "_", cmd_arguments$classB, "_", cmd_arguments$chrA, "_", cmd_arguments$chrB, "_", cmd_arguments$genomeA, "_", cmd_arguments$genomeB, ".png", collapse = "")), width = 6000, height = 6000, pointsize = 80)

    par(mar = c(4, 4, 4, 4), oma = c(1, 1, 1, 1))
    hors = hors[order(hors$SNV_per_kbp, decreasing = TRUE), ]

    plot(x = NULL, y = NULL,
        xlab = paste0(cmd_arguments$chrA, ", ", unit.name),
        ylab = paste0(cmd_arguments$chrA, ", ", unit.name),
        xaxt = "n", yaxt = "n",
        xlim = c(0, ax.len / unit.val),
        ylim = c(0, ax.len / unit.val),
        pch = 19, cex = 0.1,
        main = paste0("HORs", cmd_arguments$class, " ", cmd_arguments$chrA))

    colours_SNV <- colorRampPalette(c("green", "yellow", "red"))(length(hors$SNV_per_kbp)) [findInterval(hors$SNV_per_kbp, seq(0, SNV_per_kbp_in_red, length.out = length(hors$SNV_per_kbp)))]
    while (TRUE) {
      if (ax.len >= 10000000) {lwd_plot <- 1; break}
      if (ax.len >= 5000000) {lwd_plot <- 2; break}
      if (ax.len >= 2500000) {lwd_plot <- 3; break}
      if (ax.len >= 1250000) {lwd_plot <- 4; break}
      lwd_plot <- 5
      break
    }
    for (j in seq_len(nrow(hors))) {
      lines(x = c(hors$start.A.bp[j] / unit.val, hors$end.A.bp[j] / unit.val),
            y = c(hors$start.B.bp[j] / unit.val, hors$end.B.bp[j] / unit.val),
            pch = 19, lwd = lwd_plot, col = colours_SNV[j])
    }

    #x axis and vertical breaks
    axis(side = 1, las = 2, lwd.ticks = 4,
        at = bins.df$start.adjusted / unit.val,
        labels = format(round(bins.df$bins.starts / unit.val, 1), nsmall = 1))
    abline(v = bins.df$start.adjusted[1] / unit.val, lwd = 4, lty = 3)
    for (i in seq_len(nrow(bins.df))) {
      axis(side = 1, las = 2, lwd.ticks = 2,
          at = seq(bins.df$start.adjusted[i],  bins.df$end.adjusted[i], by = ticks.every)[-1] / unit.val,
          labels = format(round(seq(bins.df$bins.starts[i], bins.df$bins.ends[i], by = ticks.every)[-1] / unit.val, 1), nsmall = 1))
      abline(v = bins.df$end.adjusted[i] / unit.val, lwd = 4, lty = 3)
    }

    #y axis and horizontal breaks
    axis(side = 2, las = 2, lwd.ticks = 4,
        at = bins.df$start.adjusted / unit.val,
        labels = format(round(bins.df$bins.starts / unit.val, 1), nsmall = 1))
    abline(h = bins.df$start.adjusted[1] / unit.val, lwd = 4, lty = 3)
    for (i in seq_len(nrow(bins.df))) {
      axis(side = 2, las = 2, lwd.ticks = 2,
          at = seq(bins.df$start.adjusted[i],  bins.df$end.adjusted[i], by = ticks.every)[-1] / unit.val,
          labels = format(round(seq(bins.df$bins.starts[i], bins.df$bins.ends[i], by = ticks.every)[-1] / unit.val, 1), nsmall = 1))
      abline(h = bins.df$end.adjusted[i] / unit.val, lwd = 4, lty = 3)
    }
    dev.off()
  }

  
  cat("================================================================================\n")



  ### 08 / 00 Repetitiveness ============================================================================================
  cat(paste0(" 07 / 00 ", paste(rep(" ", 50), collapse = "")))
  cat(Sys.time())
  cat("\n")
  cat("################################################################################\n")
  repeats$hors_formed_count = 0
  repeats$hors_formed_tot_rep_normalised = 0
  for (i in seq_len(nrow(hors))) {
    repeats$hors_formed_count[hors$start_A[i] : (hors$start_A[i] + hors$block.size.in.units[i] - 1)] = repeats$hors_formed_count[hors$start_A[i] : (hors$start_A[i] + hors$block.size.in.units[i] - 1)] + 1
    repeats$hors_formed_count[hors$start_B[i] : (hors$start_B[i] + hors$block.size.in.units[i] - 1)] = repeats$hors_formed_count[hors$start_B[i] : (hors$start_B[i] + hors$block.size.in.units[i] - 1)] + 1
  }
  repeats$hors_formed_tot_rep_normalised = repeats$hors_formed_count / nrow(repeats)

  if (cmd_arguments$saveR == "y") {
    write.csv(x = repeats, file = file.path(cmd_arguments$output_folder, paste0("repeats_with_hors_", cmd_arguments$class, "_", cmd_arguments$classB, "_", cmd_arguments$chrA, "_", cmd_arguments$chrB, "_", cmd_arguments$genomeA, "_", cmd_arguments$genomeB, ".csv", collapse = "")), row.names = FALSE)
  }
  summary <- data.frame(genomeA = cmd_arguments$genomeA,
                        genomeB = cmd_arguments$genomeB,
                        chrA = cmd_arguments$chrA,
                        chrB = cmd_arguments$chrB,
                        repA = cmd_arguments$class,
                        repB = cmd_arguments$classB,
                        repA_total = split_after,
                        repB_total = nrow(repeats) - split_after,
                        A_repeats_found_in_B = sum(repeats$hors_formed_count[1 : split_after] > 0),
                        B_repeats_found_in_A = sum(repeats$hors_formed_count[split_after : nrow(repeats)] > 0))

  write.csv(x = summary, file = file.path(cmd_arguments$output_folder, paste0("summary_of_hors_", cmd_arguments$class, "_", cmd_arguments$classB, "_", cmd_arguments$chrA, "_", cmd_arguments$chrB, "_", cmd_arguments$genomeA, "_", cmd_arguments$genomeB, ".csv", collapse = "")), row.names = FALSE)
  


  cat("================================================================================\n")

  gc()
}