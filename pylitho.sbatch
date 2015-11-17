#!/bin/bash -l

###
### General job information (name, output file names)
###
#SBATCH -J pylitho
#SBATCH -e pylitho_e%j
#SBATCH -o pylitho_o%j

###
### Job run time (hh:mm:ss)
###
#SBATCH -t 2-00:00:00

###
### Hardware request information (number of compute cores, nodes, memory, etc.)
###

### Number of cores, nodes, etc.
#SBATCH -n 128
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=16

### a per-process (soft) memory limit in MB (1 GB = 1000)
#SBATCH --mem-per-cpu=10000

###
### Job queue
###
#SBATCH -p demo

###
### Script items for this run
###

ulimit -c unlimited
srun python ./pylitho.py
