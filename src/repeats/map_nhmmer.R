map_nhmmer = function(output_folder, i, representative, seqID, start, end, fasta_sequence, use_adist_scores = FALSE) {
   ## export the reference and array sequence ===============================
  repeat_file = file.path(output_folder, paste0("Array_", i, "_repeat.fasta"))
  array_file = file.path(output_folder, paste0("Array_", i, "_sequence.fasta"))
  seqinr::write.fasta(sequences = representative, names = "reference_repeat", file.out = repeat_file, open = "w")
  seqinr::write.fasta(sequences = fasta_sequence, names = seqID, file.out = array_file, open = "w")
  ## use nhmmer ============================================================
  nhmmer_table_output = file.path(output_folder, paste0("nhmmer_", i, "_output.txt"))
  system(paste("../dep/hmmer/nhmmer.exe", " ", 
              "--dna --tblout ", nhmmer_table_output, " ",
              repeat_file, " ",
              array_file, sep = ""), 
              intern = FALSE, wait = TRUE, show.output.on.console = FALSE, ignore.stdout = TRUE)
  ## read and parse the output into repeat table ===========================
  repeats_df <- read_and_format_nhmmer(nhmmer_table_output, seqID, i)
  if(nrow(repeats_df) != 0) {
  ## If set, change to edit distance based score ===========================
    if(use_adist_scores) {
      repeats_seq = lapply(seq_len(nrow(repeats_df)), function(X) paste0(fasta_sequence[repeats_df$start[X] : repeats_df$end[X]], collapse = "")[[1]])
      costs = list(insertions = 1, deletions = 1, substitutions = 1)
      repeats_df$score[repeats_df$strand == "+"] = adist(representative, repeats_seq[repeats_df$strand == "+"], costs)[1,]  / nchar(representative) * 100
      repeats_df$score[repeats_df$strand == "-"] = adist(rev_comp_string(representative), repeats_seq[repeats_df$strand == "-"])[1,]  / nchar(representative) * 100
    }
  ## adjust start and end coordinates ======================================
    repeats_df$start = repeats_df$start + start - 1
    repeats_df$end = repeats_df$end + start - 1
  }
  ## remove generated files ================================================
  file.remove(repeat_file)
  file.remove(array_file)
  file.remove(nhmmer_table_output)
  return(repeats_df)
}