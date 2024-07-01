write_align_read <- function(mafft_exe, temp_dir, sequences, name = "", options = "", seq_names = "", use_mafft = FALSE, remove_ali_mafft = TRUE, return_ali_mafft = TRUE) {
  if (!is.character(mafft_exe)) stop("mafft_exe is not a character string")
  if (!is.character(name)) stop("name is not a character string")
  if (!is.character(options)) stop("options are not a character string")
  if (inherits(seq_names, "character")) seq_names <- seq_along(sequences)

  if (options == "") options <- "--retree 2 --inputorder"

  if (use_mafft) {
    seqinr::write.fasta(sequences = as.list(sequences),
                        file.out = file.path(temp_dir, paste0(name, "temp.fasta")),
                        names = seq_names)
    system2(command = mafft_exe,
            args = paste(options, " ",
                         paste0("\"", file.path(temp_dir, paste0(name, "temp.fasta")), "\""),
                         " > ", paste0("\"", file.path(temp_dir, paste0(name, "temp.aligned.fasta")), "\""), sep = ""),
            wait = TRUE, stdout = TRUE, stderr = TRUE)
    file.remove(file.path(temp_dir, paste0(name, "temp.fasta")))

    alignment <- seqinr::read.alignment(file.path(temp_dir, paste0(name, "temp.aligned.fasta")), format = "FASTA", forceToLower = TRUE)
    seqinr::write.fasta(sequences = alignment$seq, 
                        file.out = file.path(temp_dir, paste0(name, "temp.aligned.fasta")),
                        names = alignment$nam)
    # Sys.sleep(0.1)
    if (remove_ali_mafft) file.remove(file.path(temp_dir, paste0(name, "temp.aligned.fasta")))
    if (return_ali_mafft) {
      return(alignment) 
    } else {
      return()
    }
  }

  sequences_bs <- Biostrings::DNAStringSet(sequences)
  names(sequences_bs) <- seq_names
  msa_message <- capture.output({alignment <- msa::msa(sequences_bs, method = "ClustalOmega", type = "dna",
                                                       verbose = FALSE)})
  # if(msa_message != "use default substitution matrix") {
  if (msa_message != "using Gonnet") {
    print(msa_message)
  }
  return(alignment)
}
