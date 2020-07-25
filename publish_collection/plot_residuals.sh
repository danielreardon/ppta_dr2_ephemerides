#!/usr/bin/tcsh

setenv dr2dir /fred/oz002/dreardon/ppta_dr2_ephemerides/publish_collection/dr2
setenv dr2edir /fred/oz002/dreardon/ppta_dr2_ephemerides/publish_collection/dr2e/ecliptic

cd $dr2edir
mkdir -p $dr2edir/output/

foreach parfile (`ls *.par | grep -v "kop_alt" | grep -v "0437"`)
	echo $parfile
	setenv psrname `echo $parfile | cut -d'.' -f1 | cut -d'_' -f1`
	tempo2 -output general2 -f $parfile ../$psrname.tim -nobs 50000 -s "{file}\t{bat}\t{freq}\t{pre}\t{post}\t{err}\t{posttn}\t{tndm}\t{tndmerr}\t{tnrn}\t{tnrnerr}\t{tnchrom}\t{tnchromerr}\n" > output/$parfile.out

end

rm temp.par

cd ..
