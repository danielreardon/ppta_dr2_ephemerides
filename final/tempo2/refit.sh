#!/usr/bin/tcsh

foreach psr ( `ls -1 *.par | cut -d'.' -f1` )
#foreach psr ( `cat psrs.list` )
	echo $psr	
	#tempo2 -f $psr.par $psr.tim -nobs 35000 -newpar
	#mv new.par $psr.new.par
	tempo2 -output general2 -f $psr.par $psr.tim -nobs 50000 -s "{file}\t{bat}\t{freq}\t{pre}\t{post}\t{err}\t{posttn}\t{tndm}\t{tndmerr}\t{tnrn}\t{tnrnerr}\t{tnchrom}\t{tnchromerr}\n" > $psr.out
end

