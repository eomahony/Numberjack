#! /bin/bash
#This script needs only the folder name as a parameter

_data="$1/*.dzn"
_model="$1/*.mzn"
echo "DATA files : "
for d in $_data
do
      echo "$d"
done


echo "Models : "
for m in $_model
do
      echo "$m"
done

output=$1"_fzn"
mkdir $output
chmod +w $output

for m in $_model
do 
for d in $_data
do
#      echo " mzn2fzn $m $d  "
l_m=`expr length $m` 
l_d=`expr length $d` 
l_1=`expr length $1` 
output_file=${m}"_"${d}".fzn" 
#echo $output_file
#echo  mzn2fzn $m $d -O- -o $output/$output_file
     mzn2fzn $m $d -O- -o $output/$output_file
echo "generated file --> $output_file"
done
done 

