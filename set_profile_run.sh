#!/bin/bash

fname="$1x$2"
ncpus=$(($1 * $2))
mem=$((4 * $ncpus))
echo $fname $ncpus $mem

echo "creating job file $fname.job"
cp 4x60.job "$fname.job"
sed -i "s/4x60/$fname/g; s/ncpus=240/ncpus=$ncpus/; s/mem=960GB/mem=${mem}GB/" "$fname.job"

#echo "creating job file ${fname}_ini.job"
#cp 4x60_ini.job "${fname}_ini.job"
#sed -i "s/4x60/$fname/g; s/ncpus=240/ncpus=$ncpus/; s/mem=960GB/mem=${mem}GB/" "${fname}_ini.job"

echo "creating in file $fname.in"
cp ROMS/External/4x60.in "ROMS/External/$fname.in"
sed -i "s/.*NtileI ==.*/NtileI == $1/; s/.*NtileJ ==.*/NtileJ == $2/" "ROMS/External/$fname.in"

#echo "creating in file ${fname}_ini.in"
#cp ROMS/External/4x60.in "ROMS/External/${fname}_ini.in"
#sed -i "s/.*NtileI ==.*/NtileI == $1/; s/.*NtileJ ==.*/NtileJ == $2/" "ROMS/External/${fname}_ini.in"



echo "submitting jobs $fname.job"
qsub $fname.job
#qsub ${fname}_ini.job

