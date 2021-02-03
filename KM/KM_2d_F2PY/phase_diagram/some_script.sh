#!/bin/bash
shopt -s expand_aliases
source ~/.aliasrc
source ~/.functionsrc
ARGUMENTS=${@}
PATH_TO_ME=${0}
SCRIPT=${0##*/}
LOGFILE=log.${SCRIPT%.sh}
>$LOGFILE
exec >  >(tee -a $LOGFILE)
exec 2> >(tee -a $LOGFILE >&2)
echo "Executing: $SCRIPT from $PATH_TO_ME"


rm -fv  *VSUloc*.dat 
while read lam U M Sx Sy Rx Ry a;do    
    echo $lam $Sx $Sy >> Sx_SyVSlam_U$U.dat
    echo $lam $Rx $Ry >> Rx_RyVSlam_U$U.dat
done <params.run


for file in *VSlam_U*;do
    sort -n $file > tmp
    mv tmp $file;
done


#plot 'params.run' u 1:2:6 with points palette pointsize 3 pointtype 5, 2-0.5*x
