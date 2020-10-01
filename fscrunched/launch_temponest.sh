#!/bin/bash
#SBATCH --job-name=tnest_J1600
#SBATCH --output=/fred/oz002/dreardon/ppta_dr2_ephemerides/fscrunched/tnest.out
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=5g
#SBATCH --tmp=4G

source ~/.cshrc_tempo2_rms

cd /fred/oz002/dreardon/ppta_dr2_ephemerides/fscrunched
srun -N1 -n1 cp -r /fred/oz002/rshannon/tempo2 $JOBFS/tempo2
export TEMPO2=$JOBFS/tempo2

export PGPLOT_DIR=/fred/oz002/rshannon/pgplot
export LD_LIBRARY_PATH=/fred/oz002/rshannon/lib:$LD_LIBRARY_PATH
export TEMPO2=/fred/oz002/rshannon/tempo2
export TEMPO2_CLOCK_DIR=/home/dreardon/tempo2/clock

srun /fred/oz002/rshannon/bin/tempo2 -gr temponest -f J1600-3053.par J1600-3053.selected.tim -cfile White-Red-DM.cfile

# White noise and Red noise and interesting model parameters
# srun tempo2 -gr temponest -f data/dr2e/J0437-4715.dr2e.par data/dr2e/J0437-4715.dr2e.tim -cfile configs/White-Band-Group-DM-Model.cfile -nobs 35000

# White	noise and Red noise
# srun tempo2 -gr temponest -f data/dr2e/J0437-4715_dr2e.par data/dr2e/J0437-4715_dr2e.tim -cfile configs/White-Band-Group-DM.cfile -nobs 35000

# Fixed white noise with Red noise and model
# srun tempo2 -gr temponest -f data/tnest_white/White-Band-Group-DM-J0437-4715-.par data/dr2e/J0437-4715_dr2e.tim -cfile configs/Band-Group-DM-Model.cfile -nobs 35000

# Fixed white noise with Red noise and model all parameters
#  srun tempo2 -gr temponest -f data/tnest_white/White-Band-Group-DM-J0437-4715-.par data/dr2e/J0437-4715_dr2e.tim -cfile configs/Band-Group-DM-Model-All.cfile -nobs 35000
