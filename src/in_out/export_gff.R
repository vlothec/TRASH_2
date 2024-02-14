export_gff = function(annotations.data.frame = "", output = ".", file.name = "gff.table",
                      seqid = ".", source = ".", type = ".", start = "0", end = "0", score = ".", strand = ".", phase = ".", attributes = ".", 
                      attribute.names = ".") 
{
  print("Export gff function")
  print("GFF function, either add numeric value (or values) for which column(s) contains data, or text string which will be universal for that gff field")
  print("attribute.names are names added to each attribute as in attributes columns selection")
  
  if(!is.data.frame(annotations.data.frame)) stop("Provide a data frame with annotations")
  
  if(nrow(annotations.data.frame) == 0) stop("The data frame provided is does contain at least one row")
  
  output = sub("[/\\\\]$", "", output)
  
  if(!is.character(seqid)) 
    if(length(seqid) == 1) seqid = annotations.data.frame[, seqid]
    else seqid = do.call(paste, c(annotations.data.frame[, seqid], sep = "_"))
  else seqid = rep(seqid, nrow(annotations.data.frame))
  
  if(!is.character(source)) 
    if(length(source) == 1) source = annotations.data.frame[, source]
  else source = do.call(paste, c(annotations.data.frame[, source], sep = "_"))
  else source = rep(source, nrow(annotations.data.frame))
  
  if(!is.character(type)) 
    if(length(type) == 1) type = annotations.data.frame[, type]
  else type = do.call(paste, c(annotations.data.frame[, type], sep = "_"))
  else type = rep(type, nrow(annotations.data.frame))
  
  if(!is.character(start)) 
    if(length(start) == 1) start = annotations.data.frame[, start]
  else start = do.call(sum, c(annotations.data.frame[, start]))
  else start = rep(start, nrow(annotations.data.frame))
  
  if(!is.character(end)) 
    if(length(end) == 1) end = annotations.data.frame[, end]
  else end = do.call(sum, c(annotations.data.frame[, end]))
  else end = rep(end, nrow(annotations.data.frame))
  
  if(!is.character(score)) 
    if(length(score) == 1) score = annotations.data.frame[, score]
  else score = do.call(sum, c(annotations.data.frame[, score]))
  else score = rep(score, nrow(annotations.data.frame))
  
  if(!is.character(strand)) 
    if(length(strand) == 1) strand = annotations.data.frame[, strand]
  else strand = do.call(paste, c(annotations.data.frame[, strand], sep = "_"))
  else strand = rep(strand, nrow(annotations.data.frame))
  
  if(!is.character(phase)) 
    if(length(phase) == 1) phase = annotations.data.frame[, phase]
  else phase = do.call(paste, c(annotations.data.frame[, phase], sep = "_"))
  else phase = rep(phase, nrow(annotations.data.frame))
  
  if(!is.character(attributes)) 
    if(length(attributes) == 1) attributes = paste(attribute.names, annotations.data.frame[, attributes], sep = "")
  else attributes = apply(annotations.data.frame[, attributes], 1, function(x) {    paste0(attribute.names, x, collapse = ";")  })
  else attributes = rep(attributes, nrow(annotations.data.frame))
  
  
  gff_format = data.frame(seqid, # sequence ID, like chromosme name
                          source,# annotation source, like TRASH
                          type,  # annotation type, like "gene", "satellite_DNA"
                          start,
                          end,
                          score,
                          strand,# +, - or "." I think
                          phase, # 1 2 or 3 I think
                          attributes) # other things

  
  options(scipen=10)
  write.table(x = gff_format, file = paste(output, "/", file.name, ".gff", sep = ""), quote = FALSE, sep = "\t", eol = "\r", row.names = FALSE, col.names = FALSE)
  options(scipen=0)
  
}
