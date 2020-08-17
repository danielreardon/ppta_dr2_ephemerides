#!/usr/bin/tcsh


foreach psr (`ls -1 J1909*.tim | cut -d'.' -f1`)
	echo $psr
	
	foreach parfile (`ls -1 J1909-3744_ecliptic.par`)
		setenv parname `echo $parfile`

		echo $parname
		tempo2 -f $parname $psr.tim -nobs 35000 -set NITS 3 -newpar
		mv new.par $parname
	end


end

