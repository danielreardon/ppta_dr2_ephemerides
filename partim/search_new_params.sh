#!/usr/bin/tcsh

setenv datadir /fred/oz002/dreardon/ppta_dr2_ephemerides/partim/dr2e_cpsr2_prepend/

mkdir -p $datadir/new_params
cd $datadir

ls -1 J*.par | cut -d. -f1 | sort | uniq | grep -v "psrs" | grep -v "J0437-4715" | grep -v "J2241-5236" > psrs.list # do J0437 and J2241 manually

foreach psr (`cat psrs.list`)
	echo $psr
	# Get chi-squared of initial fit
	setenv chisq `tempo2 -f $psr.par $psr.tim | grep "Fit Chisq" | awk '{print $7}' | cut -d'/' -f1`
	echo "Initial chisq =  $chisq"

	cp $psr.par new.par

	while (1)
		rm -f ${psr}_fit_results.txt
		touch ${psr}_fit_results.txt

		foreach param (`cat ../params.txt`)

			if ( { grep -q $param new.par } ) then
				echo "$param in par file already"
			else
				echo $param
        			setenv chisq_new `tempo2 -f new.par $psr.tim -fit $param | grep "Fit Chisq" | awk '{print $7}' | cut -d'/' -f1`
				echo "$param $chisq_new" >> ${psr}_fit_results.txt 
			endif

		end
        
        	if ( { grep -q -E "ECC" $psr.par  } ) then
			foreach param (`cat ../t2_params.txt`)

           	     		if ( { grep -q $param new.par } ) then
                	        	echo "$param in par file already"
                		else
                       	 		echo $param
                        		setenv chisq_new `tempo2 -f new.par $psr.tim -fit $param | grep "Fit Chisq" | awk '{print $7}' | cut -d'/' -f1`
                        		echo "$param $chisq_new" >> ${psr}_fit_results.txt
                		endif

        		end
		endif

		if ( { grep -q -E "EPS1" $psr.par  } ) then
                	foreach param (`cat ../ell1_params.txt`)

                       		if ( { grep -q $param new.par } ) then
                                	echo "$param in par file already"
                        	else
                                	echo $param
                                	setenv chisq_new `tempo2 -f new.par $psr.tim -fit $param | grep "Fit Chisq" | awk '{print $7}' | cut -d'/' -f1`
                                	echo "$param $chisq_new" >> ${psr}_fit_results.txt
                        	endif

                	end
        	endif

		setenv new_param `sort -k 2 ${psr}_fit_results.txt | head -1 | awk '{print $1}'`
		setenv chisq_new `sort -k 2 ${psr}_fit_results.txt | head -1 | awk '{print $2}'`

		if (`echo "$chisq - $chisq_new > 2.0" | bc`) then
			echo "Adding $new_param"
			tempo2 -f new.par $psr.tim -fit $new_param -newpar
			grep -i "TNChromAmp" $psr.par >> new.par
			grep -i "TNSubtractChrom" $psr.par >> new.par
			setenv chisq $chisq_new
        	else
			break
		endif

	end

	cp new.par new_params/$psr.par

end
cd ..
