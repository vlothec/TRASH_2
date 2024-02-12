
load.libraries() #TODO function

run.settings = parse.settings() #TODO function

if(is.na(run.settings)) stop("Parsing the arguments failed")
if(run.settings == 1) stop()


env.kmer = 20
env.region.fraction = 0.5
env.region.threshold = 99


# run.settings is a data frame with each operation to be done in order specified in the data frame:

# 1  Name of a fasta file to run
# 2  Dependancy ID (HOR run will depend on a repeat run, if that repeat run is planned, so it will have the same Dependancy ID and will be done after the repeats)
# 3  Mode (repeats, repeats.with.families, HORs, HORs.double) others to add: synteny?, 
# 4  HOR class name (for HORs and HORs.double modes)
# 5  HOR class sequence (for HORs and HORs.double modes)
# 5  Output folder directory
# 6  Name of sequence A within the fasta file to compare (for HORs.double mode only)
# 7  Name of sequence B within the fasta file to compare (for HORs.double mode only)
# 8  max repeat size
# 9  some setting about repeat identification sensitivity (region threshold, kmer)
# 10 some setting about HOR identification sensitivity (minimum length, maximum SNP)
# 11 some setting about splitting repeats
# 12 file to write the output into, defferent for each dependancy ID
# 13 Directory of fasta file to run

initialise.parallel() #TODO function


for(dep.ID in unique(run.settings$dependancy.ID))
{
  print(paste0("Starting ID group ", dep.ID))
  run.settings.temp = run.settings[run.settings$dependancy.ID == dep.ID,]
  
  for(i in 1 : nrow(run.settings.temp))
  {
    if(run.settings.temp$mode[i] == "repeats") fn.TRASH.repeats(run.settings.temp[i]) #TODO function
    
    if(run.settings.temp$mode[i] == "repeats.with.families") fn.TRASH.repeats(run.settings.temp[i]) #TODO function
    
    if(run.settings.temp$mode[i] == "HORs") fn.TRASH.HORs(run.settings.temp[i]) #TODO function
    
    if(run.settings.temp$mode[i] == "HORs.double") fn.TRASH.HORs.double(run.settings.temp[i]) #TODO function
  }
  print(paste0("Finished ID group ", dep.ID))
}






#############
# functions #
#############

fn.TRASH.repeats = function(run.settings)
{
  sequence.and.names = read.fasta.and.list(file = file.path(run.settings[13], run.settings[1]))
  
  foreach(i = 1 : length(sequence.and.names)) %dopar% {  
    
    repetitive.regions = find.repetitive.regions(fasta = sequence.and.names[i], settings = run.settings) #TODO function
    
    if(is.na(repetitive.regions)) next()
    if(repetitive.regions == 1) next()
    if(nrow(repetitive.regions) == 0) next()
    
    write.csv(x = repetitive.regions, file = file.path(run.settings[12], "temp_out", paste0("Regions_", i, ".csv")), row.names = F)
    
    repeats = find.repeats(fasta = sequence.and.names[i], settings = run.settings, regions = repetitive.regions) #TODO function
    
    if(is.na(repeats)) next()
    if(repeats == 1) next()
    if(nrow(repeats) == 0) next()
    
    write.csv(x = repeats, file = file.path(run.settings[12], "temp_out", paste0("Repeats_", i, ".csv")), row.names = F)
    
    gc()
  }
  
  regions.all = data.frame() #TODO initialise
  repeats.all = data.frame() #TODO initialise
  
  for(i in 1 : length(sequence.and.names)) 
  {  
    if(file.exists(file.path(run.settings[12], "temp_out", paste0("Regions_", i, ".csv"))))
    {
      regions.all = rbind(regions.all, read.csv(file = file.path(run.settings[12], "temp_out", paste0("Regions_", i, ".csv"))))
    }  
    if(file.exists(file.path(run.settings[12], "temp_out", paste0("Repeats_", i, ".csv"))))
    {
      repeats.all = rbind(repeats.all, read.csv(file = file.path(run.settings[12], "temp_out", paste0("Repeats_", i, ".csv"))))
    }
  }
  
  if(nrow(repetitive.regions) != 0)
  {
    write.csv(x = repetitive.regions, file = file.path(run.settings[12], paste0("Regions_", run.settings[1], ".csv")), row.names = F)
  }
  if(nrow(repeats.all) != 0)
  {
    write.csv(x = repeats.all, file = file.path(run.settings[12], "temp_out", paste0("Repeats_", run.settings[1], ".csv")), row.names = F)
    
    plot.repeats(repeats.all, repetitive.regions) #TODO function
  }
  
  
}

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


# repetitive.regions = find.repetitive.regions(fasta = sequence.and.names[i], settings = run.settings)
find.repetitive.regions = function(fasta, settings)
{
  # 1. identify regions with low score
  
  
  
  
  
  
}


# repeats = find.repeats(fasta = sequence.and.names[i], settings = run.settings, regions = repetitive.regions)
find.repeats = function(fasta, settings, regions)
{
  
}






