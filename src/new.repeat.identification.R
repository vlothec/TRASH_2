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

k = 7 #7 is optimal, 8 slows down significantly, 9 doubles the time, 10 triples
kmers.to.check = 50 # how many top hit kmers should be checked calculating N
if(4^k < kmers.to.check) kmers.to.check = 4^k # unlikely to be a problem, but why not 
max.repeat.size = 500
min.distance.region.fraction.to.consider = 1 # [%], if none is at least that, then the top one is picked

setwd(dir = "C:/Users/vlothec/Documents/GitHub/Mikoshi-random-scripts/trash_assemblies")

fasta_full = read.fasta.and.list(file = "./da_ScuGale1.1_SUPER_15_extr2.fasta") # common truncation, insertion of other elements

fasta.region = fasta_full[[1]][13625:50798]

#################
# new method
#################

time1 = Sys.time()
counts.kmers = kcount(fasta.region, k = k)

kmer.names.default = colnames(counts.kmers)

kmer.names.new = kmer.names.default[order(counts.kmers, decreasing = T)]
counts.kmers = counts.kmers[order(counts.kmers, decreasing = T)]

distances = list()
for(i in 1 : kmers.to.check)
{
  hits = str_locate_all(string = paste(fasta.region, collapse = ""), pattern = fixed(kmer.names.new[i]))
  if(nrow(hits[[1]]) > 1)
  {
    distances = c(distances, list(hits[[1]][,1][2 : nrow(hits[[1]])] - hits[[1]][,1][1 : (nrow(hits[[1]]) - 1)])) #start next - start this
    
    print(as.numeric(names(sort(table(distances[[i]]),decreasing=TRUE)[1])))
  }
}


as.numeric(names(sort(table(unlist(distances)), decreasing=TRUE)[1:14]))
# [1] 124  94  30 248 154 218 337 475  50 351
# 124 = 94 + 30
# 248 = 124 + 124
# 154 = 124 + 30
# 218 = 124 + 94
# 337 = ?
# 475 = ?

N.table = as.data.frame(sort(table(unlist(distances)), decreasing=TRUE))
N.table$fraction = 100 * N.table$Freq/sum(N.table$Freq)

potential.N = as.numeric(as.character((N.table$Var1[N.table$fraction > 1])))
if(length(potential.N) == 0) potential.N = as.numeric(as.character((N.table$Var1[1])))

networks.N = findOperations(potential.N) # more than one network indicates multiple repeat units present

# example potential Ns: 124  94  30 248 154 218 337 475
# example networks: 
# 248 218 154 124  94  30
# 154
# 337
# 475
# in reality, it comes from a region consisting of mainly 124 bp long repeats that contain 30 bp duplication, that's why it's a most
# common distance for highest count k-mers too.
# Old versions without duplication are also present, those are the 94 bp repeats. 94 + 30 bp duplication results in 124 bp
# 154, 218 and 248 are just duplications
# note that 337 and 475 did not have networks found for them 

# N value should be:
# a member of the network, 
# as big as possible, 
# not a simple composite (so when division by all other numbers should not result in an integer)
N.values = vector(mode = "numeric", length = length(networks.N))

for(i in 1 : length(N.values))
{
  ID = 1
  N.values[i] = networks.N[[i]][ID]
  while(0 %in% (N.values[i] %% networks.N[[i]][-ID]))
  {
    ID = ID + 1
    N.values[i] = networks.N[[i]][ID]
  }
}

#################
# old method
#################
seqB = paste(fasta.region, collapse = "")

time1 = Sys.time()
kmer = k
mask.small.repeats = 4
distance = vector(mode = "numeric", length =  nchar(seqB))
for(ii in 1 : (nchar(seqB) - kmer + 1))
{
  kmer.pattern = str_sub(seqB, ii, (ii + kmer - 1))
  window.string = str_sub(seqB, (ii + 1 + mask.small.repeats), (ii + max.repeat.size))
  distance[ii] = str_locate(string = window.string, pattern = kmer.pattern)[[1]] + mask.small.repeats   # TODO change str_locate to str_match to see if faster
}
Sys.time() - time1


distance = distance[!is.na(distance)]
distance = distance[distance > mask.small.repeats]
if(length(distance) > 0)
{
  periodicity = NULL
  periodicity = new.distance.N(plot = F, distance, N.max.div = 5)
}
Sys.time() - time1


###########
# reverse subset sum problem
###########


findOperations = function(numbers) #this will check if a sume of multiple single N distances or a sum of a few of them will result in one of the values and will find unique networks of numbers with such properties 
{
  if(length(numbers) == 1) return(list(numbers))
  IDs = (1:length(numbers))[order(numbers, decreasing = F)]
  numbers = numbers[order(numbers, decreasing = F)]
  max.number = max(numbers)
  
  operations = list()
  for(i in 1 : length(numbers))
  {
    result = max.number + 1
    operation = vector(mode = "numeric", length = length(numbers)) #operation vector indicates how many each number is added to create a valid operation
    operation[i] = 1
    result = sum(rep(numbers, operation))
    #check operations composed of just one value
    while(result <= max.number)
    {
      if(result %in% numbers)
      {
        operations = c(operations, list(operation))
      } 
      operation[i] = operation[i] + 1
      result = sum(rep(numbers, operation))
    }
    #check operations composed of different values, each of them once
    n = 1
    operation = vector(mode = "numeric", length = length(numbers)) #operation vector indicates how many each number is added to create a valid operation
    operation[i] = 1
    result = sum(rep(numbers, operation))
    while(result <= max.number)
    {
      if(result %in% numbers)
      {
        operations = c(operations, list(operation))
      } 
      n = n + 1
      if(n > length(numbers)) break
      operation[n] = operation[n] + 1
      result = sum(rep(numbers, operation))
    }
  }
  #operations = operations[sapply(X = operations, FUN = sum) > 1]
  operations = unique(operations)
  
  networks = list()
  i = 0
  while(i < length(operations))
  {
    i = i + 1
    added = F
    temp.network = unique(c(sum(rep(numbers, operations[[i]])), numbers[which(operations[[i]] > 0)]))
    if(length(networks) > 0)
    {
      j = 1
      while(j <= length(networks))
      {
        if(T %in% (temp.network %in% networks[[j]]))
        {
          networks[[j]] = unique(c(networks[[j]], temp.network))
          j = length(networks)
          added = T
        }
        j = j + 1
      }
    }
    if(!added) networks = c(networks, list(temp.network))
  }
  # merge networks
  i = 1
  while(i < length(networks))
  {
    networks.to.remove = NULL
    for(j in (i+1) : length(networks))
    {
      if(T %in% (networks[[j]] %in% networks[[i]]))
      {
        networks[[i]] = unique(c(networks[[i]], networks[[j]]))
        networks.to.remove = c(networks.to.remove, j)
      }
    }
    if(!is.null(networks.to.remove)) networks = networks[-networks.to.remove]
    i = i + 1
  }
  for(i in 1 : length(networks))
  {
    networks[[i]] = networks[[i]][order(networks[[i]], decreasing = T)]
  }
  return(networks)
}




