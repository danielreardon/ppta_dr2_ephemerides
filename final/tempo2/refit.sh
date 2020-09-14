#!/usr/bin/tcsh

foreach psr ( `ls -1 *.par | cut -d'.' -f1` )
	echo $psr	
	tempo2 -f $psr.par $psr.tim -nobs 35000 -newpar
	mv new.par $psr.new.par
end

