rev_comp_string <-function(DNAstr) return(tolower(toString(Biostrings::reverseComplement(Biostrings::DNAString(DNAstr)))))
