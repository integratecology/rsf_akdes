#!/bin/bash
#SBATCH --job-name=tapir_cv        # name of the job
#SBATCH --partition=defq,intel     # partition to be used (defq, gpu or intel)
#SBATCH --time=96:00:00            # walltime (up to 96 hours)
#SBATCH --nodes=1                  # number of nodes
#SBATCH --ntasks-per-node=1        # number of tasks (i.e. parallel processes) to be started
#SBATCH --cpus-per-task=1          # number of cpus required to run the script
#SBATCH --mem-per-cpu=32G	   # memory required for process
#SBATCH --array=1-400%125    	   # set number of total simulations and number that can run simultaneously	  


module load gcc

export LD_LIBRARY_PATH="/home/alston92/software/lib64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/home/alston92/software/gdal-3.3.0/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/home/alston92/software/proj-8.0.1/lib:$LD_LIBRARY_PATH"

ldd /home/alston92/R/x86_64-pc-linux-gnu-library/3.6/terra/libs/terra.so
ldd /home/alston92/R/x86_64-pc-linux-gnu-library/3.6/rgdal/libs/rgdal.so

module load R

cd /home/alston92/proj/irsf   # where executable and data is located

date
echo "Initiating script"


if [ -f results/cv_summary_tapir.csv ]; then
	echo "Results file already exists! continuing..."
else
	echo "creating results file cv_summary.csv"
	echo "aid,ind_file,cor_irsf,cor_crsf,kld_irsf,kld_crsf,kld_r2" > results/cv_summary_tapir.csv
fi

if [ -f results/cv_data_tapir.csv ]; then
        echo "Results file already exists! continuing..."
else
        echo "creating results file cv_data.csv"
        echo "aid,prob_irsf,prob_crsf,emp_count" > results/cv_data_tapir.csv
fi


Rscript cross_validation_tapir.R ${SLURM_ARRAY_TASK_ID}     # name of script
echo "Script complete"
date
