write_align_read <- function(mafft_exe, temp_dir, sequences, name = "", options = "", seq_names = "") {
  if (!is.character(mafft_exe)) {
    stop("mafft_exe is not a character string")
  }
  if (!is.character(name)) {
    stop("name is not a character string")
  }
  if (!is.character(options)) {
    stop("options are not a character string")
  }
  if (!is.character(seq_names)) {
    stop("seq_names are not a character string")
  }

  if (seq_names == "") seq_names <- seq_along(sequences)

  if (options == "") options <- "--quiet --retree 2 --inputorder"

  seqinr::write.fasta(sequences = as.list(sequences),
                      file.out = file.path(temp_dir, paste0(name, "temp.fasta")),
                      names = seq_names)
  system(paste(mafft_exe, " ", options, " ",
               file.path(temp_dir, paste0(name, "temp.fasta")),
               " > ", file.path(temp_dir, paste0(name, "temp.aligned.fasta")), sep = ""), 
               intern = FALSE, wait = TRUE, show.output.on.console = FALSE,  ignore.stderr = TRUE, ignore.stdout = TRUE)
  alignment <- seqinr::read.alignment(file.path(temp_dir, paste0(name, "temp.aligned.fasta")), format = "FASTA", forceToLower = TRUE)
  file.remove(file.path(temp_dir, paste0(name, "temp.fasta")))
  file.remove(file.path(temp_dir, paste0(name, "temp.aligned.fasta")))
  return(alignment)
}
