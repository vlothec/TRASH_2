GetMomentVectorPS <- function(seq) {
  n <- nchar(seq)
  
  # Binary indicator sequences of A, C, G, T
  uA <- numeric(n)
  uC <- numeric(n)
  uG <- numeric(n)
  uT <- numeric(n)
  
  # Sequences' length of A, C, G, T respectively
  nA <- 0
  nC <- 0
  nG <- 0
  nT <- 0
  
  # Get the binary indicator sequences and their lengths
  for (i in 1:n) {
    nu <- toupper(substr(seq, i, i))
    switch(nu,
           "A" = { uA[i] <- 1; nA <- nA + 1 },
           "C" = { uC[i] <- 1; nC <- nC + 1 },
           "G" = { uG[i] <- 1; nG <- nG + 1 },
           "T" = { uT[i] <- 1; nT <- nT + 1 })
  }
  
  # Discrete Fourier Transforms
  UA <- fft(uA)
  UC <- fft(uC)
  UG <- fft(uG)
  UT <- fft(uT)
  
  # Exclude the first term
  UA <- UA[-1]
  UC <- UC[-1]
  UG <- UG[-1]
  UT <- UT[-1]
  
  # Get the first half of the transform (since it's symmetric)
  m <- floor((n - 1) / 2)
  UA1 <- UA[1:m]
  UC1 <- UC[1:m]
  UG1 <- UG[1:m]
  UT1 <- UT[1:m]
  
  # Power spectrums
  PSA <- abs(UA)^2
  PSC <- abs(UC)^2
  PSG <- abs(UG)^2
  PST <- abs(UT)^2
  
  # Normalized moments
  MA <- numeric(3)
  MC <- numeric(3)
  MG <- numeric(3)
  MT <- numeric(3)
  
  # Moment vector
  for (j in 1:3) {
    MA[j] <- (nA * (n - nA)) * sum(PSA^j) / (nA^j * (n - nA)^j)
    MC[j] <- (nC * (n - nC)) * sum(PSC^j) / (nC^j * (n - nC)^j)
    MG[j] <- (nG * (n - nG)) * sum(PSG^j) / (nG^j * (n - nG)^j)
    MT[j] <- (nT * (n - nT)) * sum(PST^j) / (nT^j * (n - nT)^j)
  }
  
  v <- c(MA[1], MC[1], MG[1], MT[1], MA[2], MC[2], MG[2], MT[2], MA[3], MC[3], MG[3], MT[3])
  
  return(v)
}


getEDistance <- function(A, B) {
  # Calculate squared differences
  diffSquare <- (A - B)^2
  
  # Calculate the sum of squared differences along rows and take the square root
  EDist <- sqrt(rowSums(diffSquare))
  
  return(EDist)
}



# library(ape)


seqs = seqinr::read.fasta("C:/Users/vlothec/Desktop/3 documents from hash.fasta", 
                          as.string = T, seqonly = T, forceDNAtolower = T)

# seqs = seqs[[1]]
# string_length <- nchar(seqs)
# start_substrings <- lapply(seq_len(string_length), function(X) substr(seqs, 1, X))
# end_substrings <- lapply(seq_len(string_length), function(X) substr(seqs, (X + 1), string_length))
# shifted_strings <- lapply(seq_len(string_length), function(X) paste0(end_substrings[X], start_substrings[X]))
# shifted_strings_reverse <- unlist(lapply(shifted_strings, rev_comp_string))
# shifted_strings_all = c(shifted_strings, shifted_strings_reverse)
# seqs = shifted_strings_all

len <- length(seqs)

lenX <- numeric(len)
for (i in 1:len) {
  lenX[i] <- nchar(seqs[[i]])
}
lenX

b <- max(lenX)
cat('Min:', min(lenX), '\n')
cat('Max:', max(lenX), '\n')

timeA = Sys.time()
v <- list()
for (i in 1:len) {
  # Concatenate the sequences before passing to GetMomentVectorPS
  v[[i]] <- 1/b * GetMomentVectorPS(paste(seqs[[i]], collapse = ""))
}
v_matrices <- lapply(v, function(vec) matrix(vec, nrow = 1))
# Get (Euclidean) lower triangular distance matrix based on above moment vectors
D <- matrix(0, nrow = len, ncol = len)
for (j in 1:len) {
  for (i in j:len) {
    D[i, j] <- getEDistance(v_matrices[[i]], v_matrices[[j]])
    # Populate the lower triangular part symmetrically
    D[j, i] <- D[i, j]
  }
}
# Rearrange the above distance matrix into a row vector in order to use seqlinkage
print(Sys.time() - timeA)


timeA = Sys.time()
dists = adist(seqs[1], seqs)
print(Sys.time() - timeA)


mafft_dists = c(0,18,14,17,14,20,15,18,15,12,7,18,10,22,6,20,15,8,9,7,15,11)

plot(D[1,])
plot(dists[1,])
plot(mafft_dists)
cor(D[1,], dists[1,])
cor(D[1,], mafft_dists)
cor(dists[1,], mafft_dists)








