map_nhmmer = function(output_folder, arrayID, representative, seqID, start, end, fasta_sequence, nhmmer_dir) {
   ## export the reference and array sequence ===============================
  repeat_file = file.path(output_folder, paste0("Array_", arrayID, "_repeat.fasta"))
  array_file = file.path(output_folder, paste0("Array_", arrayID, "_sequence.fasta"))
  seqinr::write.fasta(sequences = representative, names = "reference_repeat", file.out = repeat_file, open = "w")
  seqinr::write.fasta(sequences = fasta_sequence, names = seqID, file.out = array_file, open = "w")
  ## use nhmmer ============================================================
  nhmmer_table_output = file.path(output_folder, paste0("nhmmer_", arrayID, "_output.txt"))
  system2(command = nhmmer_dir, 
          args = paste("--popen 0.1 --pextend 0.8 --dna --tblout ", nhmmer_table_output, " ", repeat_file, " ", array_file, sep = ""),
          wait = TRUE, stdout = TRUE, stderr = TRUE)
  ## read and parse the output into repeat table ===========================
  repeats_df <- read_and_format_nhmmer(nhmmer_table_output, seqID, arrayID)
  if(nrow(repeats_df) != 0) {
  ## adjust start and end coordinates ======================================
    repeats_df$start = repeats_df$start + start - 1
    repeats_df$end = repeats_df$end + start - 1
  }
  # write.csv(repeats_df, file.path(output_folder, "nhmmer.csv"))
  ## remove generated files ================================================
  unlink(repeat_file)
  unlink(array_file)
  unlink(nhmmer_table_output)
  remove(representative, fasta_sequence, nhmmer_table_output)
  gc()
  return(repeats_df)
}