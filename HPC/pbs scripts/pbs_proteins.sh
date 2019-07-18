#PBS -l walltime=01:00:00

#PBS -l select=1:ncpus=1:mem=1gb



module load anaconda3/personal


cd ~/OMICS project



Rscript Proteins_script.R
