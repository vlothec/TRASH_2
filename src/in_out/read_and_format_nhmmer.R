read_and_format_nhmmer <- function(nhmmer_file, seqID, arrayID) { # nolint
  if (!file.exists(nhmmer_file)) warning(paste0("read_and_format_nhmmer Warning: could not find the nhmmer file: ",
                                                nhmmer_file))
  if (file.size(nhmmer_file) == 0) warning(paste0("read_and_format_nhmmer Warning: the nhmmer file is empty ",
                                                  nhmmer_file))
  lines <- readLines(con = nhmmer_file)
  if (length(lines) == 0) {
    return(data.frame(seqID = vector(mode = "character"),
                      arrayID = vector(mode = "numeric"),
                      start = vector(mode = "numeric"),
                      end = vector(mode = "numeric"),
                      strand = vector(mode = "character"),
                      score = vector(mode = "numeric",),
                      eval = vector(mode = "numeric")))
  }
  lines <- lines[!grepl("#", x = lines)]
  if (length(lines) == 0) {
    return(data.frame(seqID = vector(mode = "character"),
                      arrayID = vector(mode = "numeric"),
                      start = vector(mode = "numeric"),
                      end = vector(mode = "numeric"),
                      strand = vector(mode = "character"),
                      score = vector(mode = "numeric",),
                      eval = vector(mode = "numeric")))
  }

  data <- lapply(lines, function(line) {
    strsplit(trimws(line), "\\s+")[[1]]
  })

  # Convert the list to a data frame
  nhmmer <- as.data.frame(do.call(rbind, data))

  if (nrow(nhmmer) == 0) {
    return(data.frame(seqID = vector(mode = "character"),
                      arrayID = vector(mode = "numeric"),
                      start = vector(mode = "numeric"),
                      end = vector(mode = "numeric"),
                      strand = vector(mode = "character"),
                      score = vector(mode = "numeric",),
                      eval = vector(mode = "numeric")))
  }

  colnames(nhmmer) <- c("target name", "accession", "query name", "accession",
                        "hmmfrom", "hmm to", "alifrom", "ali to", "envfrom",
                        "env to", "sq len", "strand", "E-value", "score",
                        "bias", "description of target")

  backup_start <- nhmmer$alifrom
  nhmmer$alifrom[nhmmer$strand == "-"] <- nhmmer$`ali to`[nhmmer$strand == "-"]
  nhmmer$`ali to`[nhmmer$strand == "-"] <- backup_start[nhmmer$strand == "-"]

  backup_start <- nhmmer$envfrom
  nhmmer$envfrom[nhmmer$strand == "-"] <- nhmmer$`env to`[nhmmer$strand == "-"]
  nhmmer$`env to`[nhmmer$strand == "-"] <- backup_start[nhmmer$strand == "-"]

  nhmmer$seqID <- seqID
  nhmmer$arrayID <- arrayID

  nhmmer <- nhmmer[c("seqID", "arrayID", "envfrom", "env to", "strand", "score", "E-value")]
  names(nhmmer) <- c("seqID", "arrayID", "start", "end", "strand", "score", "eval")

  nhmmer$seqID <- as.character(nhmmer$seqID)
  nhmmer$start <- as.numeric(nhmmer$start)
  nhmmer$end <- as.numeric(nhmmer$end)
  nhmmer$score <- as.numeric(nhmmer$score)
  nhmmer$eval <- as.numeric(nhmmer$eval)

  return(nhmmer)

}