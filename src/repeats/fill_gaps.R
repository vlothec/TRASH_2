fill_gaps <- function(repeat_table, fasta_content, array_start) {
  gap_range <- 0.5
  max_score <- 0.8
  if (nrow(repeat_table) < 2) {
    return(repeat_table[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class", "representative", "score_template")])
  }
  mean_rep_len <- round(mean(repeat_table$width))
  one_rep_range <- floor(mean_rep_len * (1 - gap_range)) : ceiling(mean_rep_len * (1 + gap_range))
  two_rep_range <- mean_rep_len + floor(mean_rep_len * (1 - gap_range)) : ceiling(mean_rep_len * (1 + gap_range))
  i <- 2
  while (i <= nrow(repeat_table)) {
    gap_len <- repeat_table$start[i] - repeat_table$end[i - 1] - 1
    if (gap_len %in% one_rep_range) {
      new_start <- repeat_table$end[i - 1] + 1
      new_end <- repeat_table$start[i] - 1
      print(new_start)
      check_seq_fw <- paste0(fasta_content[(new_start - array_start + 1) : (new_end - array_start + 1)], collapse = "")
      check_seq_rv <- rev_comp_string(check_seq_fw)
      score <- adist(repeat_table$representative[i], c(check_seq_fw, check_seq_rv))[1, ] / nchar(repeat_table$representative[i])
      if (sum(score < max_score) > 0) {
        repeat_table[nrow(repeat_table) + 1, ] <- list(seqID = repeat_table$seqID[i], 
                                                       arrayID = repeat_table$arrayID[i],
                                                       start = new_start,
                                                       end = new_end,
                                                       strand = c("+", "-")[which.min(score)],
                                                       score = -1,
                                                       eval = -1,
                                                       width = new_end - new_start + 1,
                                                       class = repeat_table$class[i],
                                                       representative = repeat_table$representative[i],
                                                       score_template = min(score))
        repeat_table <- repeat_table[order(repeat_table$start, decreasing = FALSE), ]
        i <- i + 1
      }
    }
    i <- i + 1
  }
  return(repeat_table)
}