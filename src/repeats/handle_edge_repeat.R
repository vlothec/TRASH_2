handle_edge_repeat <- function(repeats_df, sequence_vector, sequence_vector_start, template_sequence = "") {
  # repeats within an array are checked for gaps and for each edge repeat (at least one on each edge,
  # but also possible internal gaps), the edge ones are checked for best end position, and if they
  # have a sufficiently high score, potential extra repeat is checked. This check is calculating
  # local similarity

  ### dev settings
  # sequence_vector_full <- read_fasta_and_list("C:\\Users\\Piotr Włodzimierz\\Documents\\GitHub\\TRASH_dev/testing_fastas/chr5_extraction.fasta")
  # sequence_vector <- sequence_vector_full[[1]][5051 : 200000]
  # sequence_vector_start <- 5050 + 1

  # repeats_df_full <- read.csv(file = "C:\\Users\\Piotr Włodzimierz\\Documents\\GitHub\\TRASH_dev\\temp/chr5_extraction.fasta_repeats.csv")
  # repeats_df <- repeats_df_full#[repeats_df_full$start > 13943869 & repeats_df_full$start < 15189763, ]
  # repeats_df$representative <- tolower("caagcttcttcttgcttctcaaagctttgatggtgtagtcgtagtccgtatgagtctttgtctttgtatcttctaacaaggatacaatacttaggcttttaagatccgattgcggttctagttgttatactcactcatacacatgacatctagtcatatttgactccaaaacactaac")
  # template_sequence <- repeats_df$representative[1]
  ####
  potential_seq_width <- 1.2 # 120% of repeat width will be checked to find the best edge annotation size

  if (nrow(repeats_df) < 2) {
    return(repeats_df)
  }
  sequence_vector_start <- sequence_vector_start + 1

  # print(paste(sequence_vector_start, paste(sequence_vector[1:50], collapse = ""), unique(repeats_df$representative)))

  array_representative <- strsplit(repeats_df$representative[1], split = "")[[1]]
  array_representative_rev <- strsplit(tolower(toString(Biostrings::complement(Biostrings::DNAString(paste(array_representative, collapse = ""))))), split = "")[[1]]
  mean_rep_size <- mean(repeats_df$width)

  left_edges <- repeats_df$start[1]
  left_edges_rep_id <- 1
  right_edges <- NULL
  right_edges_rep_id <- NULL
  for (i in 2 : (nrow(repeats_df))) {
    if ((repeats_df$start[i] - repeats_df$end[i - 1]) > mean_rep_size) {
      left_edges <- c(left_edges, repeats_df$start[i])
      left_edges_rep_id <- c(left_edges_rep_id, i)
      right_edges <- c(right_edges, repeats_df$end[i - 1])
      right_edges_rep_id <- c(right_edges_rep_id, (i - 1))
    }
  }
  right_edges <- c(right_edges, repeats_df$end[nrow(repeats_df)])
  right_edges_rep_id <- c(right_edges_rep_id, nrow(repeats_df))

  if (length(left_edges) != length(right_edges)) stop("handle_edge_repeat: wrong edge definition")

  ### 1. Adjust edge repeat length for a maximum alignment score with the representative

  for (i in seq_along(left_edges)) {
    # print(i)
    ## left edge repeat
    if (repeats_df$strand[left_edges_rep_id[i]] == "+") {
      # if on plus strand
      # adjust the existing one
      potential_start1 <- (repeats_df$end[left_edges_rep_id[i]] - sequence_vector_start + 1 - round(mean_rep_size * potential_seq_width, 0))
      potential_end1 <- (repeats_df$end[left_edges_rep_id[i]] - sequence_vector_start + 1)
      # print(paste(potential_start1, potential_end1))
      if(potential_start1 < 1) potential_start1 = 1
      if(potential_end1 < 1) potential_end1 = 1
      sequence_potential <- rev(sequence_vector[potential_start1 : potential_end1])
      repeats_df$start[left_edges_rep_id[i]] <- repeats_df$end[left_edges_rep_id[i]] - find_edge_best_start_end(sequence_potential, rev(array_representative))
      repeats_df$eval[left_edges_rep_id[i]] <- -1
      repeats_df$score[left_edges_rep_id[i]] <- adist(repeats_df$representative[left_edges_rep_id[i]], paste0(sequence_vector[(repeats_df$start[left_edges_rep_id[i]] - sequence_vector_start + 1) : potential_end1], collapse = ""))[1, ]  / nchar(repeats_df$representative[left_edges_rep_id[i]]) * 100
      if (template_sequence != "") repeats_df$score_template[left_edges_rep_id[i]] <- adist(template_sequence, paste0(sequence_vector[(repeats_df$start[left_edges_rep_id[i]] - sequence_vector_start + 1) : potential_end1], collapse = ""))[1, ]  / nchar(template_sequence) * 100
      left_edges[i] <- repeats_df$start[left_edges_rep_id[i]]
      # and test for the new one after it, unless the previous one is already under half length
      if((repeats_df$end[left_edges_rep_id[i]] - repeats_df$start[left_edges_rep_id[i]] + 1) >= (mean_rep_size / 2)) {
        potential_start3 <- (repeats_df$start[left_edges_rep_id[i]] - sequence_vector_start + 1 - round(mean_rep_size * potential_seq_width, 0))
        potential_end3 <- (repeats_df$start[left_edges_rep_id[i]] - sequence_vector_start)
      # print(paste(potential_start3, potential_end3))
        if(potential_start3 < 1) potential_start3 = 1
        if(potential_end3 < 1) potential_end3 = 1
        sequence_potential <- rev(sequence_vector[potential_start3 : potential_end3])
        best_new_start <- find_edge_best_start_end(sequence_potential, rev(array_representative))
        if (best_new_start > 0) {
          repeats_df <- rbind(repeats_df, repeats_df[left_edges_rep_id[i], ])
          repeats_df$start[nrow(repeats_df)] <- repeats_df$start[left_edges_rep_id[i]] - best_new_start - 1
          repeats_df$end[nrow(repeats_df)] <- repeats_df$start[left_edges_rep_id[i]] -  1
          repeats_df$score[nrow(repeats_df)] <- adist(repeats_df$representative[left_edges_rep_id[i]], paste0(sequence_vector[repeats_df$start[nrow(repeats_df)] : repeats_df$end[nrow(repeats_df)]], collapse = ""))[1, ]  / nchar(repeats_df$representative[left_edges_rep_id[i]]) * 100
          if (template_sequence != "") repeats_df$score_template[nrow(repeats_df)] <- adist(template_sequence, paste0(sequence_vector[(repeats_df$start[nrow(repeats_df)] - sequence_vector_start + 1) : (repeats_df$end[nrow(repeats_df)] - sequence_vector_start + 1)], collapse = ""))[1, ]  / nchar(template_sequence) * 100
          left_edges[i] <- repeats_df$start[nrow(repeats_df)]
        }
      }
    } else {
      # if on minus strand
      # adjust the existing one
      potential_start1 <- (repeats_df$end[left_edges_rep_id[i]] - sequence_vector_start + 1 - round(mean_rep_size * potential_seq_width, 0))
      potential_end1 <- (repeats_df$end[left_edges_rep_id[i]] - sequence_vector_start + 1)
      # print(paste(potential_start1, potential_end1))
      if(potential_start1 < 1) potential_start1 = 1
      if(potential_end1 < 1) potential_end1 = 1
      sequence_potential <- rev(sequence_vector[potential_start1 : potential_end1])
      repeats_df$start[left_edges_rep_id[i]] <- repeats_df$end[left_edges_rep_id[i]] - find_edge_best_start_end(sequence_potential, array_representative_rev)
      repeats_df$eval[left_edges_rep_id[i]] <- -1
      repeats_df$score[left_edges_rep_id[i]] <- adist(repeats_df$representative[left_edges_rep_id[i]], paste0(sequence_vector[(repeats_df$start[left_edges_rep_id[i]] - sequence_vector_start + 1) : potential_end1], collapse = ""))[1, ]  / nchar(repeats_df$representative[left_edges_rep_id[i]]) * 100
      if (template_sequence != "") repeats_df$score_template[left_edges_rep_id[i]] <- adist(template_sequence, paste0(sequence_vector[(repeats_df$start[left_edges_rep_id[i]] - sequence_vector_start + 1) : potential_end1], collapse = ""))[1, ]  / nchar(template_sequence) * 100
      left_edges[i] <- repeats_df$start[left_edges_rep_id[i]]
      # and test for the new one after it, unless the previous one is already under half length
      if((repeats_df$end[left_edges_rep_id[i]] - repeats_df$start[left_edges_rep_id[i]] + 1) >= (mean_rep_size / 2)) {
        potential_start3 <- (repeats_df$start[left_edges_rep_id[i]] - sequence_vector_start + 1 - round(mean_rep_size * potential_seq_width, 0))
        potential_end3 <- (repeats_df$start[left_edges_rep_id[i]] - sequence_vector_start)
      # print(paste(potential_start3, potential_end3))
        if(potential_start3 < 1) potential_start3 = 1
        if(potential_end3 < 1) potential_end3 = 1
        sequence_potential <- rev(sequence_vector[potential_start3 : potential_end3])
        best_new_start <- find_edge_best_start_end(sequence_potential, array_representative_rev)
        if (best_new_start > 0) {
          repeats_df <- rbind(repeats_df, repeats_df[left_edges_rep_id[i], ])
          repeats_df$start[nrow(repeats_df)] <- repeats_df$start[left_edges_rep_id[i]] - best_new_start - 1
          repeats_df$end[nrow(repeats_df)] <- repeats_df$start[left_edges_rep_id[i]] -  1
          repeats_df$score[nrow(repeats_df)] <- adist(repeats_df$representative[left_edges_rep_id[i]], paste0(sequence_vector[repeats_df$start[nrow(repeats_df)] : repeats_df$end[nrow(repeats_df)]], collapse = ""))[1, ]  / nchar(repeats_df$representative[left_edges_rep_id[i]]) * 100
          if (template_sequence != "") repeats_df$score_template[nrow(repeats_df)] <- adist(template_sequence, paste0(sequence_vector[(repeats_df$start[nrow(repeats_df)] - sequence_vector_start + 1) : (repeats_df$end[nrow(repeats_df)] - sequence_vector_start + 1)], collapse = ""))[1, ]  / nchar(template_sequence) * 100
          left_edges[i] <- repeats_df$start[nrow(repeats_df)]
        }
      }
    }
    ## Right edge repeat
    if (repeats_df$strand[right_edges_rep_id[i]] == "+") {
      # if on plus strand
      potential_start2 <- (repeats_df$start[right_edges_rep_id[i]] - sequence_vector_start + 1)
      potential_end2 <- (repeats_df$start[right_edges_rep_id[i]] - sequence_vector_start + round(mean_rep_size * potential_seq_width, 0))
      # print(paste(potential_start2, potential_end2))
      if(potential_start2 > length(sequence_vector)) potential_start2 = length(sequence_vector)
      if(potential_end2 > length(sequence_vector)) potential_end2 = length(sequence_vector)
      sequence_potential <- sequence_vector[potential_start2 : potential_end2]
      repeats_df$end[right_edges_rep_id[i]] <- repeats_df$start[right_edges_rep_id[i]] + find_edge_best_start_end(sequence_potential, array_representative)
      repeats_df$eval[right_edges_rep_id[i]] <- -1
      repeats_df$score[right_edges_rep_id[i]] <- adist(repeats_df$representative[right_edges_rep_id[i]], paste0(sequence_vector[potential_start2 : (repeats_df$end[right_edges_rep_id[i]] - sequence_vector_start + 1)], collapse = ""))[1, ]  / nchar(repeats_df$representative[right_edges_rep_id[i]]) * 100
      if (template_sequence != "") repeats_df$score_template[right_edges_rep_id[i]] <- adist(template_sequence, paste0(sequence_vector[potential_start2 : (repeats_df$end[right_edges_rep_id[i]] - sequence_vector_start + 1)], collapse = ""))[1, ]  / nchar(template_sequence) * 100
      right_edges[i] <- repeats_df$end[right_edges_rep_id[i]]
      # and test for the new one after it, unless the previous one is already under half length
      if((repeats_df$end[right_edges_rep_id[i]] - repeats_df$start[right_edges_rep_id[i]] + 1) >= (mean_rep_size / 2)) {
        potential_start4 <- (repeats_df$end[right_edges_rep_id[i]] - sequence_vector_start + 2)
        potential_end4 <- (repeats_df$end[right_edges_rep_id[i]] - sequence_vector_start + 1 + round(mean_rep_size * potential_seq_width, 0))
      # print(paste(potential_start4, potential_end4))
        if(potential_start4 > length(sequence_vector)) potential_start4 = length(sequence_vector)
        if(potential_end4 > length(sequence_vector)) potential_end4 = length(sequence_vector)
        sequence_potential <- sequence_vector[potential_start4 : potential_end4]
        best_new_end <- find_edge_best_start_end(sequence_potential, array_representative)
        if (best_new_end > 0) { # only if the new start is equal to 1, is it discarded, proper filtering is in the function above
          repeats_df <- rbind(repeats_df, repeats_df[right_edges_rep_id[i], ])
          repeats_df$start[nrow(repeats_df)] <- repeats_df$end[right_edges_rep_id[i]] + 1  
          repeats_df$end[nrow(repeats_df)] <- repeats_df$end[right_edges_rep_id[i]] + best_new_end + 1
          repeats_df$score[nrow(repeats_df)] <- adist(repeats_df$representative[right_edges_rep_id[i]], paste0(sequence_vector[(repeats_df$start[nrow(repeats_df)] - sequence_vector_start + 1) : (repeats_df$end[nrow(repeats_df)] - sequence_vector_start + 1)], collapse = ""))[1, ]  / nchar(repeats_df$representative[right_edges_rep_id[i]]) * 100
          if (template_sequence != "") repeats_df$score_template[nrow(repeats_df)] <- adist(template_sequence, paste0(sequence_vector[(repeats_df$start[nrow(repeats_df)] - sequence_vector_start + 1) : (repeats_df$end[nrow(repeats_df)] - sequence_vector_start + 1)], collapse = ""))[1, ]  / nchar(template_sequence) * 100
          right_edges[i] <- repeats_df$end[nrow(repeats_df)]
        }
      }
    } else {
      # if on minus strand
      potential_start2 <- (repeats_df$start[right_edges_rep_id[i]] - sequence_vector_start + 1)
      potential_end2 <- (repeats_df$start[right_edges_rep_id[i]] - sequence_vector_start + round(mean_rep_size * potential_seq_width, 0))
      # print(paste(potential_start2, potential_end2))
      if(potential_start2 > length(sequence_vector)) potential_start2 = length(sequence_vector)
      if(potential_end2 > length(sequence_vector)) potential_end2 = length(sequence_vector)
      sequence_potential <- sequence_vector[potential_start2 : potential_end2]
      repeats_df$end[right_edges_rep_id[i]] <- repeats_df$start[right_edges_rep_id[i]] + find_edge_best_start_end(sequence_potential, rev(array_representative_rev))
      repeats_df$eval[right_edges_rep_id[i]] <- -1
      repeats_df$score[right_edges_rep_id[i]] <- adist(rev_comp_string(repeats_df$representative[right_edges_rep_id[i]]), paste0(sequence_vector[potential_start2 : (repeats_df$end[right_edges_rep_id[i]] - sequence_vector_start + 1)], collapse = ""))[1, ]  / nchar(repeats_df$representative[right_edges_rep_id[i]]) * 100
      if (template_sequence != "") repeats_df$score_template[right_edges_rep_id[i]] <- adist(rev_comp_string(template_sequence), paste0(sequence_vector[potential_start2 : (repeats_df$end[right_edges_rep_id[i]] - sequence_vector_start + 1)], collapse = ""))[1, ]  / nchar(template_sequence) * 100
      right_edges[i] <- repeats_df$end[right_edges_rep_id[i]]
      # and test for the new one after it, unless the previous one is already under half length
      if((repeats_df$end[right_edges_rep_id[i]] - repeats_df$start[right_edges_rep_id[i]] + 1) >= (mean_rep_size / 2)) {
        potential_start4 <- (repeats_df$end[right_edges_rep_id[i]] - sequence_vector_start + 2)
        potential_end4 <- (repeats_df$end[right_edges_rep_id[i]] - sequence_vector_start + 1 + round(mean_rep_size * potential_seq_width, 0))
      # print(paste(potential_start4, potential_end4))
        if(potential_start4 > length(sequence_vector)) potential_start4 = length(sequence_vector)
        if(potential_end4 > length(sequence_vector)) potential_end4 = length(sequence_vector)
        sequence_potential <- sequence_vector[potential_start4 : potential_end4]
        best_new_end <- find_edge_best_start_end(sequence_potential, rev(array_representative_rev))
        if (best_new_end > 1) { # only if the new start is equal to 1, is it discarded, proper filtering is in the function above
          repeats_df <- rbind(repeats_df, repeats_df[right_edges_rep_id[i], ])
          repeats_df$start[nrow(repeats_df)] <- repeats_df$end[right_edges_rep_id[i]] + 1  # TODO: add or subtract one here?
          repeats_df$end[nrow(repeats_df)] <- repeats_df$end[right_edges_rep_id[i]] + best_new_end + 1
          repeats_df$score[nrow(repeats_df)] <- adist(rev_comp_string(repeats_df$representative[right_edges_rep_id[i]]), paste0(sequence_vector[(repeats_df$start[nrow(repeats_df)] - sequence_vector_start + 1) : (repeats_df$end[nrow(repeats_df)] - sequence_vector_start + 1)], collapse = ""))[1, ]  / nchar(repeats_df$representative[right_edges_rep_id[i]]) * 100
          if (template_sequence != "") repeats_df$score_template[nrow(repeats_df)] <- adist(rev_comp_string(template_sequence), paste0(sequence_vector[(repeats_df$start[nrow(repeats_df)] - sequence_vector_start + 1) : (repeats_df$end[nrow(repeats_df)] - sequence_vector_start + 1)], collapse = ""))[1, ]  / nchar(template_sequence) * 100
          right_edges[i] <- repeats_df$end[nrow(repeats_df)]
        }
      }
    }
  }
  # print("handled")

  repeats_df$width <- repeats_df$end - repeats_df$start + 1
  repeats_df <- repeats_df[repeats_df$width > 2,]
  # print(nrow(repeats_df))

  ### dev settings
  # write.csv(repeats_df, "C:/Users/Piotr Włodzimierz\\Documents\\GitHub\\CC-Col\\v2_5s_adjusting/added_repeats_5s.csv")


  # repeats_df$start <- repeats_df$start - sequence_vector_start + 1
  # repeats_df$end <- repeats_df$end - sequence_vector_start + 1
  # export_gff(annotations.data.frame = repeats_df, output = "C:/Users/Piotr Włodzimierz/Documents/GitHub/TRASH_dev/temp", file.name = "added_repeats_5s.gff",
  #                     seqid = "chr5_extr_2", source = "TRASH", type = "adjusted_satellites", start = 3, end = 4, score = ".", strand = 5, phase = ".", attributes = ".", 
  #                     attribute.names = ".")




  repeats_df <- repeats_df[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class", "representative", "score_template")]
}