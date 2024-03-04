map_default = function(i, representative, seqID, start, fasta_sequence) {

  repeats_df = NULL
  max_mismatch = ceiling(nchar(representative) / 10)
  ## find forward
  match_fw <- Biostrings::matchPattern(pattern = representative, subject = fasta_sequence, max.mismatch = max_mismatch)
  if (length(match_fw) > 0) {
    match_fw = as.data.frame(match_fw)
    match_fw$start[match_fw$start < 1] = 1
    match_fw$end[match_fw$end > nchar(fasta_sequence)] = nchar(fasta_sequence)
    match_fw$strand = "+"
    match_fw$seqID = seqID
    match_fw$arrayID = i
    match_fw$start = match_fw$start + start - 1
    match_fw$end = match_fw$end + start - 1
    match_fw$score = 0
    match_fw$eval = -1
    repeats_df = rbind(repeats_df, match_fw)
  }

  ## find reverse
  match_rev <- Biostrings::matchPattern(pattern = rev_comp_string(representative), subject = fasta_sequence, max.mismatch = max_mismatch)
  if (length(match_rev) > 0) {
    match_rev = as.data.frame(match_rev)
    match_rev$start[match_rev$start < 1] = 1
    match_rev$end[match_rev$end > nchar(fasta_sequence)] = nchar(fasta_sequence)
    match_rev$strand = "-"
    match_rev$seqID = seqID
    match_rev$arrayID = i
    match_rev$start = match_rev$start + start - 1
    match_rev$end = match_rev$end + start - 1
    match_rev$score = 0
    match_rev$eval = -1
    repeats_df = rbind(repeats_df, match_rev)
  }

  if(inherits(repeats_df, "data.frame")) {
    repeats_df <- repeats_df[c("seqID", "arrayID", "start", "end", "strand", "score", "eval")]
    return(repeats_df)
  } 
  return(data.frame(seqID = vector(mode = "character"),
                    arrayID = vector(mode = "numeric"),
                    start = vector(mode = "numeric"),
                    end = vector(mode = "numeric"),
                    strand = vector(mode = "character"),
                    score = vector(mode = "numeric",),
                    eval = vector(mode = "numeric")))
  
  
}