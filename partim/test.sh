#!/usr/bin/tcsh

setenv dr1e_dir /fred/oz002/dreardon/ppta_dr2_ephemerides/partim/dr1e
setenv dr2_dir /fred/oz002/dreardon/ppta_dr2_ephemerides/partim/dr2_boris
setenv outdir /fred/oz002/dreardon/ppta_dr2_ephemerides/partim/dr2e_cpsr2_replace 

cd $dr1e_dir
ls -1 | cut -d. -f1 | sort | uniq | grep -v "psrs" | grep -v "J0437-4715" > psrs.list # do J0437 manually
cd ../..

setenv outdir /fred/oz002/dreardon/ephemeris_ppta_dr2/partim/dr2e_cpsr2_prepend

echo "PREPENDING DR1 CPSR2_20CM data, with dr1"
foreach psr (`cat $dr1e_dir/psrs.list`)

        echo $psr
        cp $dr2_dir/$psr.par $outdir/.
        cp $dr2_dir/$psr.tim $outdir/.

	setenv first_cpsr2 `grep -i "cpsr2" $dr2_dir/$psr.tim | grep "20cm" | grep -v "C " | awk '{print $3}' | uniq | sort | head -1`

        awk -v mjd="$first_cpsr2" '$3 < mjd' $dr1e_dir/$psr.tim | grep -i "cpsr2" | grep -i "20cm\|multi\|h-oh" | grep -v "T2E" > $outdir/extra.tim
        awk -v mjd="$first_cpsr2" '$3 < mjd' $dr1e_dir/$psr.tim | grep -i "cpsr2" $dr1e_dir/$psr.par | grep -i "20cm\|multi\|h-oh"  > $outdir/extra.par
        grep -i "cpsr1" $dr1e_dir/$psr.tim | grep -v "T2E" >> $outdir/extra.tim
        grep -i "cpsr1" $dr1e_dir/$psr.par >> $outdir/extra.par
        grep -i "fptm" $dr1e_dir/$psr.tim | grep -v "T2E" >> $outdir/extra.tim
        grep -i "fptm" $dr1e_dir/$psr.par >> $outdir/extra.par
        grep -i "_s21" $dr1e_dir/$psr.tim | grep -v "T2E" >> $outdir/extra.tim
        grep -i "_s21" $dr1e_dir/$psr.par >> $outdir/extra.par
        grep -i "_s22" $dr1e_dir/$psr.tim | grep -v "T2E" >> $outdir/extra.tim
        grep -i "_s22" $dr1e_dir/$psr.par >> $outdir/extra.par

        sed -i 's/ -g / -x dr1 -group /g' $outdir/extra.tim
        sed -i 's/ 7 -f / pks -f /g' $outdir/extra.tim
        sed -i 's/ -g / -group /g' $outdir/extra.par

        grep -v "B " $outdir/extra.tim > legacy.tim
        grep "B " $outdir/extra.tim > withB.tim
        rm $outdir/extra.tim
        sed -i 's/ -x dr1 / -x dr1 -B 20CM /g' legacy.tim
        cat withB.tim >> legacy.tim
        mv legacy.tim $outdir/extra.tim
        rm withB.tim

	cd $outdir
	tempo2 -f $psr.par $psr.tim -nofit -fit F0 -fit F1 -fit DM -fit FD1 -newpar
	cd ..

        mkdir -p $outdir/early_only
    	grep -v "TN" $outdir/new.par | grep -v "START" | grep -v "FINISH" > $outdir/early_only/$psr.par
        touch $outdir/early_only/$psr.tim
        echo "FORMAT 1"    >> $outdir/early_only/$psr.tim
        echo "MODE 1" >> $outdir/early_only/$psr.tim
        cat $outdir/extra.tim >> $outdir/early_only/$psr.tim

        cat $outdir/extra.tim >> $outdir/$psr.tim
        cat $outdir/extra.par >> $outdir/$psr.par

end

rm $outdir/extra.tim
rm $outdir/extra.par