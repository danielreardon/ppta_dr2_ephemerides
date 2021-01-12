#!/bin/bash
#SBATCH --job-name=J0437
#SBATCH --output=/fred/oz002/dreardon/ppta_dr2_ephemerides/partim/dr2_boris_20200810/J0437.out
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=20g
#SBATCH --tmp=4G

source ~/.cshrc_tempo2_rms

cd /fred/oz002/dreardon/ppta_dr2_ephemerides/partim/dr2_boris_20200810
srun -N1 -n1 cp -r /fred/oz002/rshannon/tempo2 $JOBFS/tempo2
export TEMPO2=$JOBFS/tempo2

export PGPLOT_DIR=/fred/oz002/rshannon/pgplot
export LD_LIBRARY_PATH=/fred/oz002/rshannon/lib:$LD_LIBRARY_PATH
export TEMPO2=/fred/oz002/rshannon/tempo2
export TEMPO2_CLOCK_DIR=/home/dreardon/tempo2/clock

/fred/oz002/rshannon/bin/tempo2 -f J0437-4715.par J0437-4715.tim -set NITS 2 -newpar -nobs 40000

