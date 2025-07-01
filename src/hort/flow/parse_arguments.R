parse_arguments <- function(arguments, run_dir) {

  spec <- matrix(c("output_folder", "o", 1, "character",
                   "hor_threshold", "t", 2, "integer", # def is 10, as in 10% of divergence
                   "hor_min_len", "l", 2, "integer",   # def is 3, as in at least 3 repeats to form a HOR
                   "class", "c", 1, "character",
                   "repeats", "r", 1, "character",
                   "method", "m", 1, "integer",       # TODO: sanitise this
                   "chrA", "A", 1, "character",
                   "chrB", "B", 2, "character",       # TODO: sanitise this
                   "repeatsB", "b", 2, "character",   # TODO: sanitise this
                   "classB", "C", 2, "character",   # TODO: sanitise this
                   "genomeA", "g", 1, "character",   # TODO: sanitise this
                   "genomeB", "G", 2, "character",   # TODO: sanitise this
                   "saveR", "s", 2, "character",   # TODO: sanitise this
                   "plot_simple", "p", 2, "character"),    # TODO: sanitise this
                 ncol = 4, byrow = TRUE)

  arg_options <- getopt::getopt(spec)

  if (is.null(arg_options$hor_threshold)) {
    arg_options$hor_threshold <- 25
  }
  if (is.null(arg_options$hor_min_len)) {
    arg_options$hor_min_len <- 3
  }
  if (is.null(arg_options$class)) {
    stop("No class provided, use -c")
  }
  if (is.null(arg_options$output_folder)) {
    print("No output folder file provided, defaulting to this directory")
    arg_options$output_folder <- run_dir
  }
  if (!dir.exists(arg_options$output_folder)) {
    stop("Specified output folder does not exist")
  }
  if (is.null(arg_options$repeats)) {
    stop("No repeat file provided")
  }
  if (!file.exists(arg_options$repeats)) {
    stop("Specified repeats file does not exist")
  }
  if (is.null(arg_options$chrA)) {
    stop("Specify the sequence name to be analysed with -A")
  }
  return(arg_options)
}