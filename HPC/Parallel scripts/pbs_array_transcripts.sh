#!/bin/sh


#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=24:mem=10gb
#PBS -N pbs_array_parallel_transcripts




cd ~/OMICSproject
module load anaconda3/personal
nchunks=5






Rscript array_parallel_transcripts_script_HPC.R $nchunks
