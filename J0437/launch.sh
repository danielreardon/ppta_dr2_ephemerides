#!/bin/bash
#SBATCH --job-name=J0437_dr2
#SBATCH --output=/fred/oz002/dreardon/ppta_dr2_ephemerides/J0437/job_dr2.out
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=20g

export PGPLOT_DIR=/fred/oz002/rshannon/pgplot
export LD_LIBRARY_PATH=/fred/oz002/rshannon/lib:$LD_LIBRARY_PATH
export TEMPO2=/fred/oz002/rshannon/tempo2
export TEMPO2_CLOCK_DIR=/home/dreardon/tempo2/clock

cd /fred/oz002/dreardon/ppta_dr2_ephemerides/J0437
/fred/oz002/rshannon/bin/tempo2 -f new.dr2.par J0437-4715.tim -nobs 35000 -newpar
mv new.par new.dr2.par
#/fred/oz002/rshannon/bin/tempo2 -output general2 -f J0437-4715.dr2e.par J0437-4715.dr2e.tim -nobs 35000 -s "{file}\t{bat}\t{freq}\t{pre}\t{post}\t{err}\t{posttn}\t{tndm}\t{tndmerr}\t{tnrn}\t{tnrnerr}\t{tnchrom}\t{tnchromerr}\t{binphase}\n" > J0437-4715.dr2e.output
