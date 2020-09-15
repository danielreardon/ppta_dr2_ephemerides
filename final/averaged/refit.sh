#!/usr/bin/tcsh

foreach psr ( `ls -1 *.par | cut -d'.' -f1` )
#foreach psr ( `cat psrs.list` )
	echo $psr	
	tempo2 -f $psr.par $psr.avg.tim -newpar
	mv new.par $psr.new.par
end

