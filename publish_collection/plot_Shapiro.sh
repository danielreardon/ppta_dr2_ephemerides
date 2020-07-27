#!/usr/bin/tcsh

setenv dr2dir /fred/oz002/dreardon/ppta_dr2_ephemerides/publish_collection/dr2/ecliptic/
setenv dr2edir /fred/oz002/dreardon/ppta_dr2_ephemerides/publish_collection/dr2e/

cd $dr2dir
mkdir -p $dr2dir/shapiro/

#foreach parfile (`grep -l "SINI \|STIG " *.par | grep -v "kop_alt" | grep -v "1545\|2145"`)
foreach parfile (`grep -l "SINI \|STIG " *.par | grep "1125"`)
	echo $parfile
	setenv psrname `echo $parfile | cut -d'.' -f1`
	tempo2 -output general2 -f $parfile $psrname.tim -nobs 50000 -s "{file}\t{bat}\t{freq}\t{pre}\t{post}\t{err}\t{posttn}\t{tndm}\t{tndmerr}\t{tnrn}\t{tnrnerr}\t{tnchrom}\t{tnchromerr}\t{binphase}\n" > shapiro/$parfile.out
	grep -v "H3 \|H4 \|STIG \|M2 " $parfile > temp.par
	tempo2 -output general2 -f temp.par $psrname.tim -nobs 50000 -s "{file}\t{bat}\t{freq}\t{pre}\t{post}\t{err}\t{posttn}\t{tndm}\t{tndmerr}\t{tnrn}\t{tnrnerr}\t{tnchrom}\t{tnchromerr}\t{binphase}\n" > shapiro/$parfile.no_shapiro.out

end

rm temp.par

cd ..
