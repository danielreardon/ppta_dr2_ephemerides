#!/usr/bin/tcsh

setenv basedir `pwd`

setenv psr $1

setenv parfile $basedir/dr2_boris/new_params_ver1/$psr.par
setenv timfile $basedir/dr2_boris/$psr.tim
setenv outfile $basedir/dr2_boris/$psr.kopeikin.out 

rm -f $outfile
touch $outfile
echo "KOM KIN KOM_FIT KIN_FIT CHISQ" > $outfile

grep -v "XDOT" $parfile  > $basedir/$psr.test.par
setenv parfile $basedir/$psr.test.par

foreach kom (`seq 0 10 350`)
	foreach kin (`seq -1 0.05 1 | awk '{print atan2(sqrt(1-$1*$1), $1)*180/3.1415926535}'`)
		echo "testing KOM=$kom, KIN=$kin"
		tempo2 -f $parfile $timfile -set KIN $kin -set KOM $kom -fit KIN -fit KOM -set NITS 3 > $psr.test
		setenv chisq `cat $psr.test | grep "Fit Chisq" | awk '{print $7}' | tail -1 | cut -d'/' -f1`
		setenv kom_fit `cat $psr.test | grep " Y" | grep "KOM" | awk '{print $3}' | tail -1`
        	setenv kin_fit `cat $psr.test | grep " Y" | grep "KIN" | awk '{print $3}' | tail -1`
		echo "$kom $kin $kom_fit $kin_fit $chisq" >> $outfile
end
