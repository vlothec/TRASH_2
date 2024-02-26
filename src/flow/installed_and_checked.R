installed_and_checked <- function() {
  suppress_messages <- TRUE
  if (!dir.exists("../R_libs")) {
    dir.create("../R_libs")
  }
  .libPaths(c(.libPaths(), gsub("src", "R_libs", getwd())))

  if (!requireNamespace("stringr", quietly = TRUE)) {
    install.packages("stringr", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  }
  if (!requireNamespace("stringdist", quietly = TRUE)) {
    install.packages("stringdist", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  }
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    install.packages("seqinr", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    # install.packages("doParallel", lib = "../R_libs", repos = "http://cran.us.r-project.org")
    install.packages("doParallel", repos = "http://cran.us.r-project.org") #This needs to be in the main directory for the parallel workers to have access to it I believe # TODO: check if true
  }
  # if (!requireNamespace("doSNOW", quietly = TRUE)) {
  #   install.packages("doSNOW", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  # }
  if (!requireNamespace("getopt", quietly = TRUE)) {
    install.packages("getopt", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    install.packages("ape", lib = "../R_libs", repos = "http://cran.us.r-project.org")
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    install.packages("BiocManager", lib = "../R_libs", repos = "http://cran.us.r-project.org")
    BiocManager::install("Biostrings")
  }
  if (suppress_messages) {
    suppressMessages(library(stringr))
    suppressMessages(library(stringdist))
    suppressMessages(library(seqinr))
    suppressMessages(library(doParallel))
    suppressMessages(library(getopt))
    suppressMessages(library(ape))
    suppressMessages(library(Biostrings))
  } else {
    library(stringr)
    library(stringdist)
    library(seqinr)
    library(doParallel)
    library(getopt)
    library(ape)
    library(Biostrings)
  }
  return(TRUE)
}
# on Windows, to compile and use HMMER, you need to install Cygwin (with what packages), compile HMMER and then use via (source? cygwin terminal?)
# 1. create a makefile file according to the manual 
# 2. Edit the file changing the prefix to somewhere it can exist (is it needed?)
# 3. Run make install
# nhmmscan   /cygdrive/c/cygwin64/home/vlothec/hmmer