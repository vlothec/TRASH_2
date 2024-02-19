read_and_format_nhmmer = function(nhmmer_file, seqID, arrayID) {

  lines <- readLines(con = nhmmer_file)
  lines <- lines[!grepl("#", x = lines)]
  if(length(lines) == 0) {
    return(data.frame(seqID = vector(mode = "numeric"), 
                      arrayID = vector(mode = "numeric"), 
                      start = vector(mode = "numeric"), 
                      end = vector(mode = "numeric"), 
                      strand = vector(mode = "character")))
  }

  data <- lapply(lines, function(line) {
    strsplit(trimws(line), "\\s+")[[1]]
  })

  # Convert the list to a data frame
  nhmmer <- as.data.frame(do.call(rbind, data))


  colnames(nhmmer) <- c("target name","accession","query name","accession","hmmfrom","hmm to","alifrom","ali to",
                        "envfrom", "env to", "sq len", "strand", "E-value", "score", "bias", "description of target")

  backup_start = nhmmer$alifrom
  nhmmer$alifrom[nhmmer$strand == "-"] = nhmmer$`ali to`[nhmmer$strand == "-"]
  nhmmer$`ali to`[nhmmer$strand == "-"] = backup_start[nhmmer$strand == "-"]

  backup_start = nhmmer$envfrom
  nhmmer$envfrom[nhmmer$strand == "-"] = nhmmer$`env to`[nhmmer$strand == "-"]
  nhmmer$`env to`[nhmmer$strand == "-"] = backup_start[nhmmer$strand == "-"]

  nhmmer$seqID = seqID
  nhmmer$arrayID = arrayID

  nhmmer = nhmmer[c("seqID", "arrayID", "envfrom", "env to", "strand")]
  names(nhmmer) <- c("seqID", "arrayID", "start", "end", "strand")

  return(nhmmer)

}