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
2. Running ```TRASH.R``` for the first time will install the required R packages. See below for the required run settings

```TRASH.R``` needs to be called directly from its directory, or added to the PATH variable for easy access

If ```TRASH.R``` does not execute, add permissions by ```chmod +x ./TRASH.R```. Using ```Rscript ./TRASH.R``` might be necessary if R code is not being recognised

## Run

TRASH is run through the ```TRASH.R``` script, with fasta file and output directory arguments:

### Required run settings:
```
-o --output output directory
-f fasta file to process
```

### Optional run settings:
```
**-p --cores_no** number of cores for parallel run, default: 1
**-m --max_rep_size** maximum repeat size, default: 1000
**-i --min_rep_size** minimum repeat size, default: 7
**-t --templates** fasta file with repeat templates and their names 
```




