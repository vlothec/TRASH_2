#!/bin/bash

sbatch <<EOT
#!/bin/bash

#SBATCH -vvv
#SBATCH --partition=long 
#SBATCH --job-name=$1
#SBATCH -o /home/pwlodzimierz/ToL/Repeats_HOR_TRASH/slurm_out/$1.out.out
#SBATCH -e /home/pwlodzimierz/ToL/Repeats_HOR_TRASH/slurm_out/$1.err.err
#SBATCH --cpus-per-task=$2
#SBATCH --exclude=mikoshi-c1,mikoshi-c2,mikoshi-a1
#SBATCH --mem=40G

# exacutalbe:
source /home/pwlodzimierz/miniconda3/etc/profile.d/conda.sh
conda activate TRASH_dev
cd /home/pwlodzimierz/TRASH_dev/src
./TRASH.R -o "$4" -f "$3" -p $2
conda deactivate
echo "done"
