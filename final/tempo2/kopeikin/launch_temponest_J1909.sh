#!/bin/bash
#SBATCH --job-name=tnest_J1909
#SBATCH --output=/fred/oz002/dreardon/ppta_dr2_ephemerides/final/tempo2/kopeikin/tnest_J1909.out
#SBATCH --ntasks=8
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=5g
#SBATCH --tmp=4G

source ~/.cshrc_tempo2_rms

cd /fred/oz002/dreardon/ppta_dr2_ephemerides/final/tempo2/kopeikin
srun -N1 -n1 cp -r /fred/oz002/rshannon/tempo2 $JOBFS/tempo2
export TEMPO2=$JOBFS/tempo2

export PGPLOT_DIR=/fred/oz002/rshannon/pgplot
export LD_LIBRARY_PATH=/fred/oz002/rshannon/lib:$LD_LIBRARY_PATH
export TEMPO2=/fred/oz002/rshannon/tempo2
export TEMPO2_CLOCK_DIR=/home/dreardon/tempo2/clock

# cfile fits M2, KOM, KIN, and low-freq noise parameters
srun /fred/oz002/rshannon/bin/tempo2 -gr temponest -f J1909-3744.par J1909-3744.tim -cfile J1909-3744.cfile

