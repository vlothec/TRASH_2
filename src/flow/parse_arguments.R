parse_arguments <- function(arguments, run_dir) {
  #run settings

  spec <- matrix(c(
    "fasta_file", "f", 1, "character", # nolint
    "output_folder", "o", 1, "character",
    "cores_no", "p", 2, "integer",
    "max_rep_size", "m", 2, "integer",
    "min_rep_size", "i", 2, "integer",
    "templates", "t", 2, "character",
    "max_alignment_length", "l", 2, "integer",
    "HOR_setting_C", "c", 2, "integer",
    "HOR_setting_V", "v", 2, "integer",
    "N_max_div", "d", 2, "integer",
    "max_N_split", "n", 2, "integer",
    "smooth_percent", "s", 2, "integer"), 
  ncol = 4, byrow = TRUE)

  arg_options <- getopt::getopt(spec)

  if (is.null(arg_options$fasta_file)) {
    stop("No fasta file provided, use -f")
  }
  if (!file.exists(arg_options$fasta_file)) {
    if (!file.exists(file.path(run_dir, arg_options$fasta_file))) {
      stop("Provided fasta file not found")
    } else {
      arg_options$fasta_file <- file.path(run_dir, arg_options$fasta_file)
    }
  }
  if (is.null(arg_options$output_folder)) {
    print("No output folder file provided, defaulting to this directory")
    arg_options$output_folder <- run_dir
  }
  if (!is.null(arg_options$min_rep_size)) {
    if(arg_options$min_rep_size < 6) {
      print("TRASH is currently not handling repeats shorter than 7 bp very well,")
      print("(Especially monomic repeats will not be identified)")
      print("To facilitate better identification of short repeats,")
      print("Adjust the -m option to a smaller value")
    }
  }
  if (!dir.exists(arg_options$output_folder)) {
    stop("Specified output folder does not exist")
  }
  if (!is.null(arg_options$templates)) {
    if(!file.exists(arg_options$templates)) {
      stop("Specified template file does not exist")
    }
  }

  if (is.null(arg_options$max_rep_size)) arg_options$max_rep_size <- 1000
  if (is.null(arg_options$min_rep_size)) arg_options$min_rep_size <- 7
  if (is.null(arg_options$templates)) arg_options$templates <- 0
  if (is.null(arg_options$cores_no)) arg_options$cores_no <- 1
  if (is.null(arg_options$max_alignment_length)) arg_options$max_alignment_length <- 200000 #This can be increased, but there's no reason to align more than 200k bp of repeats # nolint
  if (is.null(arg_options$HOR_setting_C)) arg_options$HOR_setting_C <- 2
  if (is.null(arg_options$HOR_setting_V)) arg_options$HOR_setting_V <- 5
  if (is.null(arg_options$N_max_div)) arg_options$N_max_div <- 100
  if (is.null(arg_options$max_N_split)) arg_options$max_N_split <- 12
  if (is.null(arg_options$smooth_percent)) arg_options$smooth_percent <- 2

  if(!is.na(detectCores()) & ((detectCores() - 1) < arg_options$cores_no)) {
    print("TRASH warning: set amount of cores higher than the available detected, adjusting to max cores - 1")
    arg_options$cores_no <- (detectCores() - 1)
  }
  if(arg_options$cores_no < 1) arg_options$cores_no <- 1

  return(arg_options)
}