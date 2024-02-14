
calc_GC_in_a_window = function(start, sequence, window.size) {
  return(sum(sequence[start:(start+window.size-1)] %in% c("g", "c", "G", "C")) / window.size / 100)
}

calculate_GC_in_windows = function(windows.starts, sequence, bin.size = 0)
{
  if(bin.size == 0) stop("calculate.GC.in.windows: bin.size is required and it cannot be 0")
  if(length(sequence) < windows.starts[length(windows.starts)]) stop("calculate.GC.in.windows: Something's off with he sequence, either to short or not a char vector")
  return(100 * unlist(lapply(windows.starts, FUN = calc_GC_in_a_window, sequence, bin.size)))
}
