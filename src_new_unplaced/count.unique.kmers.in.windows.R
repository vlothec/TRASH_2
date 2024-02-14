library(seqinr)
library(Biostrings)
library(stringr)
library(kmer)
library(ape)

read.fasta.and.list = function(file = "")
{
  if(!file.exists(file)) {warning(paste0("File ", file, " does not exist")); return()}
  fasta_full = read.dna(file = file, format = "fasta", as.character = T) # this ape function is faster than read.fasta()
  if(typeof(fasta_full) == "character") 
  {
    names.save = rownames(fasta_full)
    fasta_full = list(fasta_full)
    names(fasta_full) = names.save
    remove(names.save)
  }
  return(fasta_full)
}

extract.kmers = function(X, kmer, fasta_sequence) return(paste(fasta_sequence[X:(X+kmer-1)], collapse = ""))

kmers.no.to.cover.P.fraction.of.sequence = function(sequence.vec = NULL, kmer = 12, fraction.P = 0.5, window.size = 10000, removeN = T)
{
  sequence.full.length = length(sequence.vec)
  starts = genomic.bins.starts(start = 1, end = sequence.full.length, bin.size = window.size)
  ends = c((starts[2:length(starts)] - 1), sequence.full.length)
  if(length(ends) != length(starts)) ends = sequence.full.length
  
  scores = vector(mode = "numeric", length = length(starts))
  for(i in 1 : length(starts))
  {
    print(paste0("Window:", i, "/", length(starts)))
    #counts.kmers = kcount(sequence.vec[starts[i]:ends[i]], k = kmer) #using kmer package, better for long windows and short kmers
    counts.kmers = table(unlist(lapply(X = starts[i]:ends[i], FUN = extract.kmers, kmer, sequence.vec))) #using own script, better for short windows and long kmers 
    
    counts.kmers = counts.kmers[!grepl("n", names(counts.kmers))]
    
    counts.kmers = counts.kmers[!grepl("N", names(counts.kmers))]
    
    if(length(counts.kmers) < kmer)
    {
      scores[i] = 100
      next
    }
    
    total.kmers = sum(counts.kmers)
  
    counts.kmers = sort(counts.kmers, decreasing = T)
    
    scores[i] = 100 * (min(which(unlist(lapply( 1:length(counts.kmers), 
                                               FUN = function(x, vector) return(sum(vector[1:x])), 
                                               counts.kmers )) / total.kmers > fraction.P))-1) / (total.kmers) / fraction.P
  }
  
  return(scores)
}


setwd(dir = "C:/Users/vlothec/Documents/GitHub/Mikoshi-random-scripts/trash_assemblies")
#setwd(dir = "C:/Users/wlodz/Documents/GitHub/Mikoshi-random-scripts/trash_assemblies")


######
# Alternative 1: Find N using kmer distances
######

sequence = "Rbrevi_chr1_extr1_edited"
{
  fasta_full = read.fasta.and.list(file = paste0("./", sequence, ".fasta")) # common truncation, insertion of other elements
  
  sequence.check = fasta_full[[1]]
  
  kmer = 20
  max.repeat = 2000
  fraction = 0.5
  threshold = 99
  
  window.size = max.repeat
  
  scores = kmers.no.to.cover.P.fraction.of.sequence(sequence.vec = sequence.check, 
                                                    kmer = kmer, 
                                                    fraction.P = fraction, 
                                                    window.size = window.size)

  #these starts and ends are calculated in the same way in the kmers.no.to.cover.P.fraction.of.sequence function, so make sure it stays that way!
  starts = genomic.bins.starts(start = 1, end = length(sequence.check), bin.size = window.size)
  ends = c((starts[2:length(starts)] - 1), length(sequence.check))
  if(length(ends) != length(starts)) ends = length(sequence.check)
  
  repetitive.regions = data.frame(starts = starts[scores < threshold], ends = ends[scores < threshold], scores = scores[scores < threshold])
  
  i = 1
  while(i < nrow(repetitive.regions))
  {
    if((repetitive.regions$ends[i] + 1) == (repetitive.regions$starts[i+1]))
    {
      repetitive.regions$scores[i] = ((repetitive.regions$scores[i] * (repetitive.regions$ends[i] - repetitive.regions$starts[i])) + 
                                        (repetitive.regions$scores[i+1] * (repetitive.regions$ends[i+1] - repetitive.regions$starts[i+1]))) / (repetitive.regions$ends[i+1] - repetitive.regions$starts[i]) 
      repetitive.regions$ends[i] = repetitive.regions$ends[i+1]
      repetitive.regions = repetitive.regions[-(i+1),]
      i = i - 1
    }
    i = i + 1
  }
  
  
  # N value calc
  
  
  kmers.to.check = 1:100 # which hit kmers should be checked calculating N
  max.repeat.size = max.repeat
  min.distance.region.fraction.to.consider = 1 # [%], if none is at least that, then the top one is picked
  
  repetitive.regions$N.k5 = 0
  repetitive.regions$consensus.seq.k5 = ""
  repetitive.regions$consensus.k5.count = 0
  repetitive.regions$N.k10 = 0
  repetitive.regions$consensus.seq.k10 = ""
  repetitive.regions$consensus.k10.count = 0
  repetitive.regions$N.k20 = 0
  repetitive.regions$consensus.seq.k20 = ""
  repetitive.regions$consensus.k20.count = 0
  repetitive.regions$N.k40 = 0
  repetitive.regions$consensus.seq.k40 = ""
  repetitive.regions$consensus.k40.count = 0
  
  for(i in 1 : nrow(repetitive.regions))
  {
    print(paste0("oeoe",i))
    
    fasta.region = fasta_full[[1]][repetitive.regions$starts[i] : repetitive.regions$ends[i]]
    string.fasta = paste(fasta.region, collapse = "")
    string.fasta.minus = tolower(revCompString(string.fasta))
    for(k in c(5,10,20,40))
    {
      col.k = which(names(repetitive.regions) == paste0("N.k", k))
      
      #counts.kmers = kcount(fasta.region, k = k)
      kmers.list = unlist(lapply(X = 1:length(fasta.region), FUN = extract.kmers, k, fasta.region))
      counts.kmers = table(kmers.list)
      
      #kmer.names.default = colnames(counts.kmers)
      kmer.names.default = names(counts.kmers)
      
      
      kmer.names.new = kmer.names.default[order(counts.kmers, decreasing = T)]
      counts.kmers = counts.kmers[order(counts.kmers, decreasing = T)]
      
      counts.kmers = counts.kmers[!grepl("N", kmer.names.new)]
      kmer.names.new = kmer.names.new[!grepl("N", kmer.names.new)]
      counts.kmers = counts.kmers[!grepl("n", kmer.names.new)]
      kmer.names.new = kmer.names.new[!grepl("n", kmer.names.new)]
      
      # find N
      
      distances = list()
      kmers.vector = NULL
      
      for(j in kmers.to.check)
      {
        hits = which(kmer.names.new[j] == kmers.list)
        #hits = str_locate_all(string = string.fasta, pattern = fixed(kmer.names.new[j])) can remove as it does not output overlapping hits
        if(length(hits) > 1)
        {
          distances = c(distances, list(hits[2:length(hits)] - hits[1:(length(hits)-1)])) #start next - start this
          kmers.vector = c(kmers.vector, kmer.names.new[j])
        }
      }
      
      if(length(distances) == 0)
      {
        repetitive.regions[i, col.k] = 0
        repetitive.regions[i, (col.k+1)] = ""
        repetitive.regions[i, (col.k+2)] = 0
        
        next
      }
      repetitive.regions[i, col.k] = as.numeric(names(sort(table(unlist(distances)), decreasing=TRUE))[1])
      
      
      if(F) #new N fixing method
      {
        if(length(unique(unlist(distances))) == 1)
        {
          repetitive.regions[i, col.k] = unique(unlist(distances))
        } else
        {
          N.table = as.data.frame(sort(table(unlist(distances)), decreasing=TRUE))
        N.table$fraction = 100 * N.table$Freq/sum(N.table$Freq)
        
        potential.N = as.numeric(as.character((N.table$Var1[N.table$fraction > 1])))
        if(length(potential.N) == 0) potential.N = as.numeric(as.character((N.table$Var1[1])))
        
        networks.N = findOperations(potential.N) 
        repetitive.regions[i, col.k] = networks.N[[1]][1]
        }
        
      }
      if(F) #old N fixing method
      {
        
      }
      if(repetitive.regions[i, col.k] < 2) 
      {
        repetitive.regions[i, col.k] = 0
        repetitive.regions[i, (col.k+1)] = ""
        repetitive.regions[i, (col.k+2)] = 0
        next
      }
        
      # find representative
      
      first.N.kmers = NULL
      for(j in 1 : length(distances))
      {
        if(repetitive.regions[i, col.k] %in% distances[[j]])
        {
          first.N.kmers = c(first.N.kmers, kmers.vector[j])
        }
      }
      
      if(length(first.N.kmers) == 0) next
      
      locations = list()
      
      match1 = NULL
      match2 = NULL
      for(j in 1 : length(first.N.kmers))
      {
        match1 = c(match1, which(first.N.kmers[j] == kmers.list))
        match2 = c(match2, which(tolower(revCompString(first.N.kmers[j])) == kmers.list))
      }
      
      match = c(match1, match2)
      
      if(length(match) == 0)
      {
        print("this should not happen, check")
        next
      }
      breaks.no = ceiling(length(fasta.region)/repetitive.regions[i, col.k])
      if(repetitive.regions[i, col.k] < (k * 2)) breaks.no = k * 2
        
      start.putative = hist(match, plot = F, 
                            breaks = breaks.no)$breaks[which.max(hist(match, plot = F, breaks = breaks.no)$counts)]
      end.putative = start.putative + repetitive.regions[i, col.k] - 1
      if(end.putative >= length(fasta.region))
      {
        start.putative = start.putative - (end.putative - length(fasta.region))
        end.putative = end.putative - (end.putative - length(fasta.region))
      }
      if(start.putative < 1) start.putative = 1
      
      repetitive.regions[i, (col.k+1)] = paste(fasta.region[start.putative:end.putative], collapse = "")
      
      #replace with hmmer
      if(repetitive.regions[i, col.k] > max.repeat.size) 
      {
        repetitive.regions[i, col.k] = 0
        repetitive.regions[i, (col.k+1)] = ""
        repetitive.regions[i, (col.k+2)] = 0
        next
      }
      repetitive.regions[i, (col.k+2)] = length(matchPattern(pattern = repetitive.regions[i, (col.k+1)], 
                                                             subject = string.fasta, 
                                                             max.mismatch = floor(repetitive.regions[i, col.k]/3))) +   
        length(matchPattern(pattern = repetitive.regions[i, (col.k+1)], 
                            subject = string.fasta.minus, 
                            max.mismatch = floor(repetitive.regions[i, col.k]/3)))
    }
    
  }
  
  repetitive.regions$k5.coverage = 100 * (repetitive.regions$consensus.k5.count * repetitive.regions$N.k5) / 
    (repetitive.regions$ends - repetitive.regions$starts)
  
  repetitive.regions$k10.coverage = 100 * (repetitive.regions$consensus.k10.count * repetitive.regions$N.k10) / 
    (repetitive.regions$ends - repetitive.regions$starts)
  
  repetitive.regions$k20.coverage = 100 * (repetitive.regions$consensus.k20.count * repetitive.regions$N.k20) / 
    (repetitive.regions$ends - repetitive.regions$starts)
  
  repetitive.regions$k40.coverage = 100 * (repetitive.regions$consensus.k40.count * repetitive.regions$N.k40) / 
    (repetitive.regions$ends - repetitive.regions$starts)
  
  repetitive.regions$top.k = 0
  repetitive.regions$top.k.N = 0
  repetitive.regions$top.k.seq = ""
  repetitive.regions$top.k.cov = 0
  repetitive.regions$top.k.count = 0
  
  for(i in 1 : nrow(repetitive.regions))
  {
    which.top.k = NULL
    which.top.k = which.max(c(repetitive.regions$k5.coverage[i], 
                              repetitive.regions$k10.coverage[i], 
                              repetitive.regions$k20.coverage[i], 
                              repetitive.regions$k40.coverage[i]))
    repetitive.regions$top.k[i] = c(5,10,20,40)[which.top.k]
    repetitive.regions$top.k.N[i] = repetitive.regions[i,grep("N.k", names(repetitive.regions))[which.top.k]]
    repetitive.regions$top.k.seq[i] = repetitive.regions[i,grep("consensus.seq.k", names(repetitive.regions))[which.top.k]]
    repetitive.regions$top.k.cov[i] = repetitive.regions[i,grep(".coverage", names(repetitive.regions))[which.top.k]]
    repetitive.regions$top.k.count[i] = repetitive.regions[i,grep(".count", names(repetitive.regions))[which.top.k]]
  }
  
  repetitive.regions.clean = repetitive.regions[, names(repetitive.regions) %in% c("starts", "ends", "scores", "top.k", "top.k.N", "top.k.seq", "top.k.cov", "top.k.count")]
  repetitive.regions.clean$ID = 1:nrow(repetitive.regions.clean)
  
  
  repeats = data.frame(ID = vector(mode = "numeric", length = 0), 
                       start = vector(mode = "numeric", length = 0), 
                       end = vector(mode = "numeric", length = 0), 
                       width = vector(mode = "numeric", length = 0), 
                       strand = vector(mode = "character", length = 0), 
                       seq.region.consensus.dir = vector(mode = "character", length = 0), 
                       region.ID = vector(mode = "numeric", length = 0))
  
  for(i in 1 : nrow(repetitive.regions.clean))
  {
    if(repetitive.regions.clean$top.k.seq[i] == "") next
    fasta.region = fasta_full[[1]][repetitive.regions.clean$starts[i] : repetitive.regions.clean$ends[i]]
    string.fasta = paste(fasta.region, collapse = "")
    
    hits.plus = matchPattern(pattern = repetitive.regions.clean$top.k.seq[i], 
                             subject = string.fasta, 
                             max.mismatch = floor(repetitive.regions.clean$top.k.N[i] / 3)) 
    if(length(hits.plus) > 0)
    {
      hits.plus = as.data.frame(hits.plus)
      hits.plus$ID = (nrow(repeats) + 1) : (nrow(repeats) + nrow(hits.plus))
      hits.plus$strand = "+"
      names(hits.plus)[grep("seq", names(hits.plus))] = "seq.region.consensus.dir"
      hits.plus$region.ID = i
      hits.plus$start = hits.plus$start + repetitive.regions.clean$starts[i]
      hits.plus$end = hits.plus$end + repetitive.regions.clean$starts[i]
      repeats = rbind(repeats, hits.plus)
      remove(hits.plus)
    }
    
    
    hits.minus = matchPattern(pattern = tolower(revCompString(repetitive.regions.clean$top.k.seq[i])), 
                              subject = string.fasta, 
                              max.mismatch = floor(repetitive.regions.clean$top.k.N[i] / 3))
    if(length(hits.minus) > 0)
    {
      hits.minus = as.data.frame(hits.minus)
      hits.minus$ID = (nrow(repeats) + 1) : (nrow(repeats) + nrow(hits.minus))
      hits.minus$strand = "-"
      names(hits.minus)[grep("seq", names(hits.minus))] = "seq.region.consensus.dir"
      hits.minus$region.ID = i
      hits.minus$start = hits.minus$start + repetitive.regions.clean$starts[i]
      hits.minus$end = hits.minus$end + repetitive.regions.clean$starts[i]
      repeats = rbind(repeats, hits.minus)
      remove(hits.minus)
    }
    
  }
  
  export.gff(annotations.data.frame = repeats, seqid = sequence, start = 1, end = 2, strand = 6, type = "satellite_DNA", file.name = sequence)
  export.gff(annotations.data.frame = repetitive.regions.clean, seqid = sequence, start = 1, end = 2, strand = ".", 
             score = 3, type = "satellite_region", file.name = paste0("regions", sequence))
  write.csv(x = repeats, file = paste0("repeats.", sequence, ".csv"), row.names = F)
  write.csv(x = repetitive.regions.clean, file = paste0("regions.clean.", sequence, ".csv"), row.names = F)
  write.csv(x = repetitive.regions, file = paste0("regions.", sequence, ".csv"), row.names = F)
}


######
# Aternative 2: Find N using kmer graphs
######

sequence = "Rbrevi_chr1_extr1_edited"
{
  fasta_full = read.fasta.and.list(file = paste0("./", sequence, ".fasta")) # common truncation, insertion of other elements
  
  sequence.check = fasta_full[[1]]
  
  kmer = 20
  max.repeat = 2000
  fraction = 0.01
  threshold = 99
  
  window.size = max.repeat + kmer + 1
  
  scores = kmers.no.to.cover.P.fraction.of.sequence(sequence.vec = sequence.check, 
                                                    kmer = kmer, 
                                                    fraction.P = fraction, 
                                                    window.size = window.size)
  #scores = 100*scores/(max.repeat*fraction)
  starts = genomic.bins.starts(start = 1, end = length(sequence.check), bin.size = window.size)
  ends = c((starts[2:length(starts)] - 1), length(sequence.check))
  
  repetitive.regions = data.frame(starts = starts[scores < threshold], ends = ends[scores < threshold], scores = scores[scores < threshold])
  
  i = 1
  while(i < nrow(repetitive.regions))
  {
    if((repetitive.regions$ends[i] + 1) == (repetitive.regions$starts[i+1]))
    {
      repetitive.regions$scores[i] = ((repetitive.regions$scores[i] * (repetitive.regions$ends[i] - repetitive.regions$starts[i])) + 
                                        (repetitive.regions$scores[i+1] * (repetitive.regions$ends[i+1] - repetitive.regions$starts[i+1]))) / (repetitive.regions$ends[i+1] - repetitive.regions$starts[i]) 
      repetitive.regions$ends[i] = repetitive.regions$ends[i+1]
      repetitive.regions = repetitive.regions[-(i+1),]
      i = i - 1
    }
    i = i + 1
  }
  
  
  # N value calc
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  repeats = data.frame(ID = vector(mode = "numeric", length = 0), 
                       start = vector(mode = "numeric", length = 0), 
                       end = vector(mode = "numeric", length = 0), 
                       width = vector(mode = "numeric", length = 0), 
                       strand = vector(mode = "character", length = 0), 
                       seq.region.consensus.dir = vector(mode = "character", length = 0), 
                       region.ID = vector(mode = "numeric", length = 0))
  
  for(i in 1 : nrow(repetitive.regions.clean))
  {
    if(repetitive.regions.clean$top.k.seq[i] == "") next
    fasta.region = fasta_full[[1]][repetitive.regions.clean$starts[i] : repetitive.regions.clean$ends[i]]
    string.fasta = paste(fasta.region, collapse = "")
    
    hits.plus = matchPattern(pattern = repetitive.regions.clean$top.k.seq[i], 
                             subject = string.fasta, 
                             max.mismatch = floor(repetitive.regions.clean$top.k.N[i] / 3)) 
    if(length(hits.plus) > 0)
    {
      hits.plus = as.data.frame(hits.plus)
      hits.plus$ID = (nrow(repeats) + 1) : (nrow(repeats) + nrow(hits.plus))
      hits.plus$strand = "+"
      names(hits.plus)[grep("seq", names(hits.plus))] = "seq.region.consensus.dir"
      hits.plus$region.ID = i
      hits.plus$start = hits.plus$start + repetitive.regions.clean$starts[i]
      hits.plus$end = hits.plus$end + repetitive.regions.clean$starts[i]
      repeats = rbind(repeats, hits.plus)
      remove(hits.plus)
    }
    
    
    hits.minus = matchPattern(pattern = tolower(revCompString(repetitive.regions.clean$top.k.seq[i])), 
                              subject = string.fasta, 
                              max.mismatch = floor(repetitive.regions.clean$top.k.N[i] / 3))
    if(length(hits.minus) > 0)
    {
      hits.minus = as.data.frame(hits.minus)
      hits.minus$ID = (nrow(repeats) + 1) : (nrow(repeats) + nrow(hits.minus))
      hits.minus$strand = "-"
      names(hits.minus)[grep("seq", names(hits.minus))] = "seq.region.consensus.dir"
      hits.minus$region.ID = i
      hits.minus$start = hits.minus$start + repetitive.regions.clean$starts[i]
      hits.minus$end = hits.minus$end + repetitive.regions.clean$starts[i]
      repeats = rbind(repeats, hits.minus)
      remove(hits.minus)
    }
    
  }
  
  export.gff(annotations.data.frame = repeats, seqid = sequence, start = 1, end = 2, strand = 6, type = "satellite_DNA", file.name = sequence)
  export.gff(annotations.data.frame = repetitive.regions.clean, seqid = sequence, start = 1, end = 2, strand = ".", 
             score = 3, type = "satellite_region", file.name = paste0("regions", sequence))
}










###############
# results of running time comparison of kmer counting using kcount and my script:
# 85 kb sequence, mostly repetitive
###############

#1kb windows
k = c(5,7,9,12) #kmer
s1 = c(1.1,7,120,1000) #average time per sequence with method 1, kcount
s2 = c(8,9,10,11) #average time per sequence with method 2, mine
plot(k,s1, ylim = c(0,15), col = "blue", type = "o")
par(new = T)
plot(k,s2, ylim = c(0,15), col = "red", type = "o")

#10kb windows
k = c(5,7,9,12)
s1 = c(0.8, 1.2, 10.5, 1000) #average time per sequence with method 1
s2 = c(8,8,9,10) #average time per sequence with method 2
plot(k,s1, ylim = c(0,15), col = "blue", type = "o")
par(new = T)
plot(k,s2, ylim = c(0,15), col = "red", type = "o")

#100kb windows (so ended up one 85 kb window)
k = c(5,7,9,12)
s1 = c(0.684401, 0.7736635, 2.08, 1000) #average time per sequence with method 1
s2 = c(9.07, 9.21, 11.01769, 9.42547) #average time per sequence with method 2
plot(k,s1, ylim = c(0,15), col = "blue", type = "o")
par(new = T)
plot(k,s2, ylim = c(0,15), col = "red", type = "o")

###############
# section end #
###############


if(F) # this is testing for counting method (kmer package and my method) and representation. My method and setting 7 are optimal
{
  
  k = 12 #7 is optimal, 8 slows down significantly, 9 doubles the time, 10 triples
  window.size = 10000
  sequence.full.length = length(fasta_full[[1]])
  
  starts = genomic.bins.starts(start = 1, end = sequence.full.length, bin.size = window.size)
  ends = c((starts[2:length(starts)] - 1), sequence.full.length)
  if(length(ends) != length(starts)) ends = sequence.full.length
  
  
  setting.text = c("setting 1 frequency of top 1 kmer", #looks good with low k and large window # looks bad with longer k and with longer windows
                   "setting 2 frequency of top Tt kmers", #looks good with low k and large window, bad with short window, maybe higher Tt helps
                   "setting 3 number of unique kmers", #kinda decent with short window
                   "setting 4 number of non-unique kmers", #decent with short window 
                   "setting 5 number of kmers with at least H hits", #adding another variable (H) is too problematic
                   "setting 6 number of kmers with less than H hits", #good with short window
                   "setting 7 min number of kmers that sum to 50% of the window") #looks good with low k and large window # best overall
  Tt = 10
  H = 10
  scores = list()
  for(i in 1 : length(setting.text)) scores = c(scores, list(NULL))
  
  timeA = Sys.time()
  for(i in 1 : length(starts))
  {
    print(paste0("Window:", i, "/", length(starts)))
    #counts.kmers = kcount(fasta_full[[1]][starts[i]:ends[i]], k = k) #using kmer package, better for long windows and short kmers
    counts.kmers = table(unlist(lapply(X = starts[i]:(ends[i]-k+1), FUN = extract.kmers, fasta_full[[1]]))) #using own script, better for short windows and long kmers 
    total.kmers = ends[i] - starts[i] - k + 2
    
    scores[[1]] = c(scores[[1]], 100 * max(counts.kmers) / total.kmers)
    scores[[2]] = c(scores[[2]], 100 * sum(sort(counts.kmers, decreasing = T)[1:Tt]) / total.kmers)
    scores[[3]] = c(scores[[3]], 100 * length(counts.kmers[counts.kmers == 1]) / total.kmers)
    scores[[4]] = c(scores[[4]], 100 * length(counts.kmers[counts.kmers > 1]) / total.kmers)
    scores[[5]] = c(scores[[5]], 100 * length(counts.kmers[counts.kmers >= H]) / total.kmers)
    scores[[6]] = c(scores[[6]], 100 * length(counts.kmers[counts.kmers < H & counts.kmers > 0]) / total.kmers)
    counts.kmers2 = sort(counts.kmers, decreasing = T)
    scores[[7]] = c(scores[[7]], min(which(unlist(lapply(2:length(counts.kmers2), 
                                                         FUN = function(x, vector) return(sum(vector[1:x])), 
                                                         counts.kmers2))/total.kmers > 0.5)))
    
  }
  timeB = Sys.time() - timeA
  for(setting.plot in 1:7)
  {
    print(setting.plot)
    png(filename = paste0("Kmer windows = ", window.size, " plot, k = ", k, " Tt = ", Tt, " H = ", H, ". ", setting.text[setting.plot], ".png"),
        width = 2000, height = 1000, pointsize = 16)
    par(mfrow = c(1,2))
    plot(NULL, xlim = c(1,sequence.full.length), ylim = c(0, max(scores[[setting.plot]])), xlab = "Coordinates, bp", ylab = "%", 
         main = paste0("W = ", window.size, " k = ", k, " Tt = ", Tt, " H = ", H, ". ", setting.text[setting.plot]))
    for(i in 1 : length(scores[[setting.plot]]))
    {
      lines(x = c(starts[i],ends[i]), y = c(scores[[setting.plot]][i],scores[[setting.plot]][i]))
    }
    hist(scores[[setting.plot]], breaks = length(scores[[setting.plot]]), main = "histogram of scores")
    dev.off()
  }
}

#plot window scores for repetitive window calc
png(filename = paste0("test.plot.png"),
    width = 2000, height = 1000, pointsize = 16)
plot(NULL, xlim = c(1,length(sequence.check)), ylim = c(0, 100), xlab = "Coordinates, bp", ylab = "%", 
     main = paste0("00"))
for(j in 1 : length(scores))
{
  lines(x = c(starts[j],ends[j]), y = c(scores[j],scores[j]))
}
dev.off()




