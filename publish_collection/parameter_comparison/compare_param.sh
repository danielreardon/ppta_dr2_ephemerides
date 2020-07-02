#!/bin/sh

for psr in `cat psr.list`; do
    echo $psr
    for param in `cat ${psr}.param`; do
        #echo $param
        awk -v param=$param '{if ($1 == param) printf "%s %s %.1e\n", $1,$2,$4}'  ../dr2/${psr}.par > dr2.val
        awk -v param=$param '{if ($1 == param) printf "%s %.1e\n", $2,$4}'  ../dr2e/${psr}.par > dr2e.val
        #grep $param ../dr2/${psr}.par | awk '{printf "%s %s %.1e\n", $1,$2,$4}' > dr2.val 
        #grep $param ../dr2e/${psr}.par | awk '{printf "%s %.1e\n",  $2,$4}' > dr2e.val 
        
         
        #paste dr2.val dr2e.val   
        paste dr2.val dr2e.val | awk '{if (1 == 1) printf "%s %s %s %s %s %.2f\n",  $1,$2,$3,$4,$5, $3/$5}' > val.dat


        echo -n "$psr "

        awk '{if (($1 != "RAJ") && ($1 != "DECJ")) printf "%s %s %s %s %s %s %.1f\n", $1,$2,$3,$4,$5,$6, ($4-$2)/sqrt($5*$5+$3*$3)}'  val.dat
        dr2val=`awk '{if (($1 == "RAJ") || ($1 == "DECJ")) print $2}' val.dat | awk -F : '{print $3}'`
        dr2eval=`awk '{if (($1 == "RAJ") || ($1 == "DECJ")) print $4}' val.dat | awk -F : '{print $3}'`  
        

        #echo $dr2val $dr2eval
        awk -v dr2val=$dr2val -v dr2eval=$dr2eval  '{if (($1 == "RAJ") || ($1 == "DECJ")) printf "%s %s %s %s %s %s %.1f\n",  $1,$2,$3,$4,$5,$6, (dr2eval-dr2val)/sqrt($5*$5+$3*$3)}' val.dat
            

    done

done
