installed_and_checked <- function() {
  suppress_messages <- TRUE
  if (!dir.exists("../R_libs")) {
    dir.create("../R_libs")
  }
  .libPaths(c(.libPaths(), gsub("src", "R_libs", getwd())))

  if (!requireNamespace("stringr", quietly = TRUE)) {
    install.packages("stringr", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  }
  # if (!requireNamespace("stringdist", quietly = TRUE)) {
  #   install.packages("stringdist", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  # }
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    install.packages("seqinr", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  }
  # if (!requireNamespace("doParallel", quietly = TRUE)) {
  #   install.packages("doParallel", repos = "http://cran.us.r-project.org") #This needs to be in the main directory for the parallel workers to have access to it I believe # TODO: check if true
  # }
  if (!requireNamespace("getopt", quietly = TRUE)) {
    install.packages("getopt", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  }
  # if (!requireNamespace("ape", quietly = TRUE)) {
  #   install.packages("ape", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  # }
  # if (!requireNamespace("Biostrings", quietly = TRUE)) {
  #   install.packages("BiocManager", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  #   BiocManager::install("Biostrings")
  # }
  # if (!requireNamespace("msa", quietly = TRUE)) {
  #   BiocManager::install("msa")
  # }
  if (suppress_messages) {
    suppressMessages(library(stringr))
    # suppressMessages(library(stringdist))
    suppressMessages(library(seqinr))
    # suppressMessages(library(doParallel))
    suppressMessages(library(getopt))
    # suppressMessages(library(ape))
    # suppressMessages(library(Biostrings))
  } else {
    library(stringr)
    # library(stringdist)
    library(seqinr)
    # library(doParallel)
    library(getopt)
    # library(ape)
    # library(Biostrings)
  }
  return(TRUE)
}