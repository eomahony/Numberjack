echo"mzn2fzn -G exec problems/fjsp/fjsp.mzn --data problems/fjsp/easy01.dzn -o alldiff.fzn --output-ozn-to-file alldiff.ozn
fzn-exec-free -a alldiff.fzn" 
mzn2fzn -G exec problems/fjsp/fjsp.mzn --data problems/fjsp/easy01.dzn -o alldiff.fzn --output-ozn-to-file alldiff.ozn
fzn-exec-free -a alldiff.fzn
grep "cumulative" alldiff.fzn
rm ./alldiff.* 
