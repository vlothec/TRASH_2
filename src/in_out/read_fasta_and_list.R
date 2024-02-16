read_fasta_and_list <- function(file = "") {
  if (!file.exists(file)) {
    stop(paste0("File ", file, " does not exist"))
  }

  fasta_full <- ape::read.dna(file = file, format = "fasta", as.character = TRUE, as.matrix = FALSE) # nolint
  if (typeof(fasta_full) == "character") {
    names.save <- rownames(fasta_full)
    for (i in seq_along(names.save)) names.save[i] <- strsplit(names.save[i], split = " ")[[1]][1]
    fasta_full <- list(fasta_full)
    names(fasta_full) <- names.save
    remove(names.save)
  }
  return(fasta_full)
}
