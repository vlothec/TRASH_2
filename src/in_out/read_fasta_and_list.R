read_fasta_and_list <- function(file = "") {
  if (!file.exists(file)) {
    stop(paste0("File ", file, " does not exist"))
  }

  fasta_full <- ape::read.dna(file = file, format = "fasta", as.character = TRUE) # nolint
  if (typeof(fasta_full) == "character") {
    names.save <- rownames(fasta_full)
    fasta_full <- list(fasta_full)
    names(fasta_full) <- names.save
    remove(names.save)
  }
  return(fasta_full)
}
