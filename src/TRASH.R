#!/usr/bin/env Rscript
# options(error = function() {traceback(2, max.lines=500); if(!interactive()) quit(save="no", status=1, runLast=T)})
timeA <- as.numeric(Sys.time())
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
source_files <- list.files(path = ".", pattern = ".R", recursive = TRUE)
source_files <- unlist(source_files[-grep("TRASH.R", source_files)])
source_files <- unlist(source_files[-grep("hort", source_files)])
source_files <- unlist(source_files[-grep("HORT.R", source_files)])
for (i in seq_along(source_files)) {
  source(source_files[i])
}

if (installed_and_checked()) {
  arguments <- parse_arguments(arguments, run_dir)
  set.seed(0)
  main(arguments)
  print("TRASH exiting correctly")
  print(paste0("Time elapsed: ", dhms(as.numeric(Sys.time()) - timeA)))
} else {
  print("TRASH exiting 1")
  print(paste0("Time elapsed: ", dhms(as.numeric(Sys.time()) - timeA)))
}
