#!/bin/bash
filecontent=( `cat "bugs.txt" `)

for t in "${filecontent[@]}"
do
echo "file == $t"
cp $t Bugs/${t:9}
done
echo "DONE"

