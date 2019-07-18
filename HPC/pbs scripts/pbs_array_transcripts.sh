#!/bin/sh


#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=24:mem=24gb
#PBS -N pbs_parallel_transcripts
#PBS -J 1-4


cd ~/OMICSproject
module load anaconda3/personal


m=$PBS_ARRAY_INDEX # model index
nchunks=24


Rscript array_parallel_transcripts_script_HPC.R $nchunks $m 
