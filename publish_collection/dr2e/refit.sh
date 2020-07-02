#!/usr/bin/tcsh


foreach psr (`ls -1 *.tim | cut -d'.' -f1`)
	echo $psr
	
	foreach parfile (`ls -1 ecliptic/$psr*par`)
		setenv parname `echo $parfile`

		echo $parname
		#echo "CLK_CORR_CHAIN pks2gps.dr2e.clk gps2utc.clk utc2tai.clk tai2tt_bipm2018.clk" >> $parname
		tempo2 -f $parname $psr.tim -nobs 30000 -set NITS 2 -newpar
		mv new.par $parname
            	echo "CLK_CORR_CHAIN pks2gps.dr2e.clk gps2utc.clk utc2tai.clk tai2tt_bipm2018.clk" >> $parname
		echo "wrote $parname"
	end


end

