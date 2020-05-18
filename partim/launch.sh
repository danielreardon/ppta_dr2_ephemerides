#!/bin/bash
#SBATCH --job-name=J2145_kopeikin
#SBATCH --output=/fred/oz002/dreardon/ppta_dr2_ephemerides/partim/J2145.out
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10g
#SBATCH --tmp=1G
cp -r $TEMPO2 $JOBFS/tempo2
cp -r $TEMPO2_CLOCK_DIR $JOBFS/tempo2_clock_dir
export TEMPO2=$JOBFS/tempo2
export TEMPO2_CLOCK_DIR=$JOBFS/tempo2_clock_dir

# Navigate to working directory
cd /fred/oz002/dreardon/ppta_dr2_ephemerides/partim/
./search_kopeikin.sh J2145-0750

