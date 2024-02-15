revCompString = function(DNAstr) 
{
  return(toupper(toString(reverseComplement(DNAString(DNAstr)))))
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

a = c(1,2,3,4,5,6)
for(i in which(a > 4)) {
  print(i)
  print(a[i])
}













