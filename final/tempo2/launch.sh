#!/bin/bash
#SBATCH --job-name=refit
#SBATCH --output=/fred/oz002/dreardon/ppta_dr2_ephemerides/final/tempo2/refit.out
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=20g

export PGPLOT_DIR=/fred/oz002/rshannon/pgplot
export LD_LIBRARY_PATH=/fred/oz002/rshannon/lib:$LD_LIBRARY_PATH
export TEMPO2=/fred/oz002/rshannon/tempo2
export TEMPO2_CLOCK_DIR=/home/dreardon/tempo2/clock

cd /fred/oz002/dreardon/ppta_dr2_ephemerides/final/tempo2
./refit.sh

