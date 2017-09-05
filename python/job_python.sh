#!/bin/bash

#SBATCH -p hernquist
#SBATCH -J sator
#SBATCH -n 1
#SBATCH -t 8640 # 6d in min
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-user=fbecerra@cfa.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#mpirun -np 32 /n/home02/fbecerra/arepo_old/sator/sator /n/home02/fbecerra/arepo_old/sator/par.txt 2 nahw1 138 138 0.00001
mpirun -np 1 python /n/home02/fbecerra/arepo_old/sator/python/test.py
#mpirun -np 16 /n/home02/fbecerra/arepo_old/sator/sator /n/home02/fbecerra/arepo_old/sator/par.txt 9 nahw1r1 138 nahw1 3 > /n/home02/fbecerra/arepo_old/sator/data/outfile_sator_$SLURM_JOB_ID.out
