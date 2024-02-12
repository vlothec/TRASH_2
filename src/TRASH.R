#!/usr/bin/env Rscript
this_file <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  match <- grep("--file=", cmd_args)
  if (length(match) > 0) {
    return(normalizePath(sub("--file=", "", cmd_args[match])))
  } else {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
arguments <- commandArgs(trailingOnly = TRUE)
run_dir <- getwd()

setwd(gsub("TRASH.R", "", this_file()))
#setwd("./src")
source_files <- list.files(path = ".", pattern = ".R", recursive = TRUE)
source_files <- unlist(source_files[-grep("TRASH.R", source_files)])
for (i in seq_along(source_files)) {
  source(source_files[i])
}

if (installed_and_checked()) {
  arguments <- parse_arguments(arguments, run_dir)
  print(arguments)
  main(arguments)
  print("Exiting 0")
} else {
  print("Exiting 1")
}