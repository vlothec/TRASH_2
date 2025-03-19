This is an early-development version of TRASH 2, update of https://github.com/vlothec/TRASH software for repeat identification

## New features:
1. Classification of repeats to repeat families/classes across the fasta file
2. Better mapping of very short and very long repeats
3. Additional polishing steps for the repeats found at array edges
4. Better parallelisation
5. Better error diagnostics and runtime progress reporting
6. Full re-write with updates to all algorithms


## Installation:
1. R needs to be installed.
2. mafft and nhmmer need to be installed.
3. Running ```TRASH.R``` from the ```/src/``` directory for the first time will install the required R packages (if they're missing). See below for the required run settings

```TRASH.R``` needs to be called directly from its directory, or added to the PATH variable for easy access

If ```TRASH.R``` does not execute, add permissions by ```chmod +x ./TRASH.R```. Using ```Rscript ./TRASH.R``` might be necessary if R code is not being recognised

mafft and nhmmer need to be installed and added to the PATH variable. Alternatively, both can be installed locally and their paths can be added to the ```src/main.R``` script, replacing lines 12 and 13 on Windows or 15 and 16 on Linux. **Windows** installation of nhmmer will require a Unix-like enviroment interface like Cygwin.


## Run

TRASH is run through the ```TRASH.R``` script founr in the ```/src/``` directory, with fasta file and output directory arguments:

### Required run settings:
```
-o --output             output directory
-f fasta                file to process
```

### Optional run settings:
```
-p --cores_no           number of cores for parallel run, default: 1
-m --max_rep_size       maximum repeat size, default: 1000
-i --min_rep_size       minimum repeat size, default: 7
-t --templates          fasta file with repeat templates and their names 
```

### Output

```bash
├── [fasta_file]
│   ├── [fasta_file]_repeats_with_seq.csv      main output file with identified repeats
│   ├── [fasta_file]_repeats.gff               main output repeat file in gff format
│   ├── [fasta_file]_repeats.csv               main output file with identified repeats without sequence column
│   ├── [fasta_file]_arrays.csv                repeat arrays, start and end are not perfectly aligned with repeats, but can be used to get locations of repeats without loading in potentially big repeat files
│   ├── [fasta_file]_arrays.gff                repeat arrays as above, in gff format
│   ├── [fasta_file]_run_time.csv              report of the script run time
│   ├── [fasta_file]_regarrays.csv             temp file, can be removed
│   ├── [fasta_file]_aregarrays.csv            temp file, can be removed
│   ├── [fasta_file]_classarrays.csv           temp file, can be removed
│   └── [fasta_file]_no_repeats_arrays.csv     temp file, can be removed
```
