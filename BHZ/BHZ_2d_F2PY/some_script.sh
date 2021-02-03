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
while read U M Tz Sz E0 Ex a;do    
    echo $U $Tz $Sz >> Tz_SzVSUloc_M$M.dat
    echo $U $E0 $Ex >> E0_ExVSUloc_M$M.dat
done <params.run


for file in *VSUloc_M*;do
    sort -n $file > tmp
    mv tmp $file;
done


#plot 'params.run' u 1:2:6 with points palette pointsize 3 pointtype 5, 2-0.5*x
