export.gff = function(annotations.data.frame = "", output = ".", file.name = "gff.table",
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


consensus_N = function(alignment, N)
{
  alignment.matrix = as.matrix.alignment(alignment)
  frequencies = vector(mode = "numeric", length = ncol(alignment.matrix))
  for(i in 1 : length(frequencies))
  {
    frequencies[i] = frequencies[i] + sum(1 * (alignment.matrix[,i] != "-"))
  }
  
  is.in.consensus = vector(mode = "logical", length = ncol(alignment.matrix))
  for(i in 1 : N)
  {
    is.in.consensus[order(frequencies, decreasing = TRUE)[i]] = TRUE
  }
  consensus = vector(mode = "logical", length = N)
  consensus.ID = 1
  for(i in 1 : ncol(alignment.matrix))
  {
    if(is.in.consensus[i])
    {
      consensus[consensus.ID] = c("g","c","t","a")[which.max(c(length(which(alignment.matrix[,i] == "g")),
                                                               length(which(alignment.matrix[,i] == "c")),
                                                               length(which(alignment.matrix[,i] == "t")),
                                                               length(which(alignment.matrix[,i] == "a"))))]
      consensus.ID = consensus.ID + 1
    }
  }
  return(toupper(paste(consensus, collapse = "")))
}


revCompString = function(DNAstr) 
{
  return(toupper(toString(reverseComplement(DNAString(DNAstr)))))
}


genomic.bins.starts = function(start = 1, end = 0, bin.number = 0, bin.size = 0)
{
  if(bin.number > 0 & bin.size > 0) stop("genomic.bins.starts: Use either bin number or bin size")
  if(bin.number == 0 & bin.size == 0) stop("genomic.bins.starts: Use either bin number or bin size")
  if(end < start) stop("genomic.bins.starts: End smaller than start will not work too well...")
  if(bin.size >= end) return(start)
  
  if(bin.number > 0)
  {
    seq.per.bin = (end - start + 1) %/% bin.number
    remaining.seq = (end - start + 1) %% bin.number
    bin.sizes = rep(seq.per.bin, bin.number)
    if(remaining.seq > 0) 
    {
      add.remaining.here = sample(1:bin.number, remaining.seq)
      bin.sizes[add.remaining.here] = bin.sizes[add.remaining.here] + 1
    }
    start.positions = bin.sizes
    for(i in 2 : length(start.positions))
    {
      start.positions[i] = start.positions[i] + start.positions[i-1]
    }
    start.positions = start.positions - start.positions[1] + start
    #if(start.positions[length(start.positions)] == end) start.positions = start.positions[-length(start.positions)]
    remove(seq.per.bin, remaining.seq)
    
    return(start.positions)
  }
  if(bin.size > 0)
  {
    start.positions = seq(start, (end - bin.size), bin.size)
    #if(bin.size * (end %/% bin.size) == end ) start.positions = start.positions[-length(start.positions)]
    if(end < bin.size) start.positions = 0
    if((end - start.positions[length(start.positions)] - bin.size) > bin.size/2) start.positions = c(start.positions, start.positions[length(start.positions)]+bin.size)
    return(start.positions)
  }
}


reorder.names.by.chr.etc = function(names.to.reorder = "")
{
  #reorder repeats so they appear in more or less correct order, start with chromosomes, then super, then atg, hap, scaf, contig and the rest
  names.to.reorder.numbers = as.numeric(gsub("\\D", "", names.to.reorder))
  which.chr = grep("chr|ch", names.to.reorder)
  which.super = grep("sup|spr", names.to.reorder)
  which.atg = grep("atg", names.to.reorder)
  which.hap = grep("hap", names.to.reorder)
  which.contig = grep("ctg|con", names.to.reorder)
  which.rest = (1:length(names.to.reorder))[!1:length(names.to.reorder) %in% c(which.chr, which.super, which.atg, which.hap, which.contig)]
  
  new.order = c(which.chr[order(names.to.reorder.numbers[which.chr])],
                which.super[order(names.to.reorder.numbers[which.super])],
                which.atg[order(names.to.reorder.numbers[which.atg])],
                which.hap[order(names.to.reorder.numbers[which.hap])],
                which.contig[order(names.to.reorder.numbers[which.contig])],
                which.rest[order(names.to.reorder.numbers[which.rest])])
  
  return(names.to.reorder[new.order])
}


calc.GC.in.a.window = function(start, sequence, window.size) {
  return(sum(sequence[start:(start+window.size-1)] %in% c("g", "c", "G", "C")) / window.size / 100)
}

calculate.GC.in.windows = function(windows.starts, sequence, bin.size = 0)
{
  if(bin.size == 0) stop("calculate.GC.in.windows: bin.size is required and it cannot be 0")
  if(length(sequence) < windows.starts[length(windows.starts)]) stop("calculate.GC.in.windows: Something's off with he sequence, either to short or not a char vector")
  return(100 * unlist(lapply(windows.starts, FUN = calc.GC.in.a.window, sequence, bin.size)))
}

calculate.repeats.percentage.in.windows = function(windows.starts, repeat.starts, repeat.lengths, sequence.length = 0)
{
  if(sequence.length < windows.starts[length(windows.starts)]) print("calculate.repeats.percentage.in.windows: sequence length is shorter than last window")
  if(sequence.length < max(repeat.starts)) print("calculate.repeats.percentage.in.windows: sequence length is shorter than last window")
  
  if(length(windows.starts) == 1) return(sum(repeat.lengths) / (sequence.length - windows.starts[1]))
  
  windows.starts.t = windows.starts[windows.starts < sequence.length]
  repeat.lengths.t = repeat.lengths[repeat.starts < sequence.length]
  repeat.starts.t = repeat.starts[repeat.starts < sequence.length]
  
  bins.breaks = c(windows.starts.t, sequence.length)
  
  hist.data = hist(rep(repeat.starts.t, repeat.lengths.t), breaks = bins.breaks, plot = F)
  
  return(100 * hist.data$counts / (bins.breaks[2:length(bins.breaks)] - bins.breaks[1:(length(bins.breaks) - 1)]))
}




create_table = function(data, headers = NULL, x_offset = 1, y_offset = 0.5, col_width = 1, row_height = 0.5, colours = "black") 
{
  # Made for visualise.families.R script in Rhynchospora analysis with OpenAI 3.5 help for coordinates calculations
  
  # Get the dimensions of the data
  num_rows = nrow(data)
  num_cols = ncol(data)
  
  colours = c("black", rep(colours, ((num_rows %/% length(colours) + 1))))
  # Increase number of rows for headers
  if(!is.null(headers)) 
  {
    data = rbind(headers, data)
    num_rows = num_rows + 1
  }
  
  # Create a blank plot with appropriate dimensions
  plot(0, 0, type = "n", xlab = "", ylab = "", xlim = c(0.5, (num_cols * col_width + 0.5)), ylim = c(0, num_rows * row_height), 
       frame.plot = F, xaxt = "n", yaxt = "n")
  
  # Add vertical lines for columns
  for (i in 1 : (num_cols-1)) abline(v = i * col_width + x_offset - col_width / 2)
  abline(h = 0)
  # Add horizontal lines for rows
  for (i in 0 : num_rows) abline(h = i * row_height + y_offset)
  
  # Add text to cells
  for (i in 1:num_rows) 
  {
    for (j in 1:num_cols) 
    {
      text(j * col_width - col_width / 2 + x_offset - col_width / 2, 
           (num_rows - i) * row_height - row_height / 2 + y_offset, 
           labels = as.character(data[i, j]), cex = 5, col = colours[i])
    }
  }
}
















