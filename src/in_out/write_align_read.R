write_align_read = function(mafft_exe, temp_dir, sequences, name = "", options = "", seq_names = "") {
  if(seq_names == "") seq_names = seq_along(sequences)
  
  if(options == "") options = "--quiet --retree 2 --inputorder"
  
  write.fasta(sequences = as.list(sequences), 
              file.out = file.path(temp_dir, paste0(name, "temp.fasta")), 
              names = seq_names)
  system(paste(mafft_exe, " ", options, " ", 
               file.path(temp_dir, paste0(name, "temp.fasta")), 
               " > ", file.path(temp_dir, paste0(name, "temp.aligned.fasta")), sep = ""), intern = TRUE)
  alignment = read.alignment(file.path(temp_dir, paste0(name, "temp.aligned.fasta")), format = "FASTA", forceToLower = TRUE)
  return(alignment)
}
