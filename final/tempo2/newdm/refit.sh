#!/usr/bin/tcsh


foreach psr (`cat psrs.list`)
	echo $psr	
	tempo2 -f $psr.par ../$psr.tim -nobs 35000 -set NITS 2 -newpar
	mv new.par $psr.new.par
end

