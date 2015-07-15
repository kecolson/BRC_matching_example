#!/bin/bash

# Job name: 

#SBATCH --job-name=match_balance
#

# Partition:

#SBATCH --partition=savio

#

# Account:

#SBATCH --account=co_biostat

#

# QoS:

#SBATCH --qos=biostat_normal
#

# Processors:

#SBATCH --ntasks=60
#

# Memory requirement:

###SBATCH --mem-per-cpu=6G

#

# Wall clock limit:

#SBATCH --time=96:00:00

#

# Mail type:

#SBATCH --mail-type=all

#

# Mail user:

#SBATCH --mail-user=kecolson@berkeley.edu

#

## Run command

module unload intel

module load gcc openmpi r Rmpi

mpirun -n 1 R --vanilla < 02_match_balance_master_with_coverage.R > 02_match_balance_master_with_coverage.Rout
