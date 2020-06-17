#!/usr/bin/tcsh

setenv dr2dir /fred/oz002/dreardon/ppta_dr2_ephemerides/publish_collection/dr2
setenv dr2edir /fred/oz002/dreardon/ppta_dr2_ephemerides/publish_collection/dr2e/

cd $dr2edir

foreach psr (`ls -1 *.tim | cut -d'.' -f1`)
	echo $psr
	
	foreach parfile (`ls -1 $dr2dir/$psr.*par`)
		setenv parname `echo $parfile | cut -d'/' -f8`
		grep -v "JUMP" $parfile | grep -v "TN" | grep -v "START" | grep -v "FINISH" > params.par
		cat $psr.par | grep "JUMP\|TN" >> params.par

		tempo2 -f params.par $psr.tim -nobs 30000 -set NITS 3 -newpar
		mv new.par $dr2edir/$parname
		echo "wrote $parname"
	end


end

rm params.par

cd ..
