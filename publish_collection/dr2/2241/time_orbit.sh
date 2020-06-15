#!/usr/bin/tcsh

setenv basedir `pwd`

setenv psr "J2241-5236"

setenv parfile $basedir/$psr.par
setenv timfile $basedir/$psr.tim
setenv outfile $basedir/$psr.tasc.txt

rm -f $outfile
touch $outfile
echo "MJD TASC TASC_ERR" > $outfile

foreach mjd (`seq 55234 50 58184`)

	setenv endmjd `echo $mjd | awk '{print $1 + 50}'`
	rm -f cut.tim
	touch cut.tim
	echo "FORMAT 1" >> cut.tim
	echo "MODE 1" >> cut.tim
	awk -v mjd="$mjd" -v endmjd="$endmjd" '$3 <= endmjd && $3 >= mjd' $psr.tim >> cut.tim
	tempo2 -f J2241-5236.par cut.tim -nofit -fit tasc  > temp.txt
	setenv tasc `cat temp.txt | grep "TASC" | awk '{print $4}'`
	setenv tascerr `cat temp.txt | grep "TASC" | awk '{print $5}'`
	setenv midmjd `echo $mjd | awk '{print $1 + 25}'`
	echo "$midmjd $tasc $tascerr"
	echo "$midmjd $tasc $tascerr" >> $outfile


end


