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

  if (options == "") options <- "--retree 2 --inputorder"

  sequences_bs <- Biostrings::DNAStringSet(sequences)
  names(sequences_bs) <- seq_names
  msa_message <- capture.output({alignment <- msa::msa(sequences_bs, method = "ClustalOmega", type = "dna",
                                                       verbose = FALSE)})
  # if(msa_message != "use default substitution matrix") {
  if (msa_message != "using Gonnet") {
    print(msa_message)
  }
  return(alignment)

  # alignment <- as.matrix(Biostrings::unmasked(alignment))
  # alignment_list = NULL
  # for(i in 1 : nrow(alignment)) {
  #   alignment_list <- append(alignment_list, list(tolower(paste(alignment[i, ], collapse = "")[[1]])))
  # }
  # return(alignment_list)


  # seqinr::write.fasta(sequences = as.list(sequences),
  #                     file.out = file.path(temp_dir, paste0(name, "temp.fasta")),
  #                     names = seq_names)
  # system2(command = mafft_exe, 
  #         args = paste(options, " ",
  #                      file.path(temp_dir, paste0(name, "temp.fasta")),
  #                      " > ", file.path(temp_dir, paste0(name, "temp.aligned.fasta")), sep = ""),
  #         wait = TRUE, stdout = TRUE, stderr = TRUE)
  # Sys.sleep(0.1)
  # alignment <- seqinr::read.alignment(file.path(temp_dir, paste0(name, "temp.aligned.fasta")), format = "FASTA", forceToLower = TRUE)
  # file.remove(file.path(temp_dir, paste0(name, "temp.fasta")))
  # file.remove(file.path(temp_dir, paste0(name, "temp.aligned.fasta")))
  # return(alignment)
}
