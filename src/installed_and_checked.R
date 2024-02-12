installed_and_checked <- function() {
  if (!dir.exists("../R_libs")) {
    dir.create("../R_libs")
  }
  .libPaths(c(.libPaths(), gsub("src", "R_libs", getwd())))

  if (!require("stringr", quietly = TRUE)) {
    install.packages("stringr", lib = "../R_libs")
  }
  if (!requireNamespace("stringdist", quietly = TRUE)) {
    install.packages("stringdist", lib = "../R_libs")
  }
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    install.packages("seqinr", lib = "../R_libs")
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    install.packages("doParallel", lib = "../R_libs")
  }
  if (!requireNamespace("getopt", quietly = TRUE)) {
    install.packages("getopt", lib = "../R_libs")
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    install.packages("ape", lib = "../R_libs")
  }
  if (!require("Biostrings", quietly = TRUE)) {
    install.packages("BiocManager", lib = "../R_libs")
    BiocManager::install("Biostrings")
  }
  library(stringr, quietly = TRUE, verbose = FALSE)
  library(stringdist, quietly = TRUE, verbose = FALSE)
  library(seqinr, quietly = TRUE, verbose = FALSE)
  library(doParallel, quietly = TRUE, verbose = FALSE)
  library(getopt, quietly = TRUE, verbose = FALSE)
  library(ape, quietly = TRUE, verbose = FALSE)
  library(Biostrings, quietly = TRUE, verbose = FALSE)

  return(TRUE)
}