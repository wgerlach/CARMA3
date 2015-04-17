#!/bin/bash
if [ $# -eq 0 ]; then
echo usage: checkChunks NUM SUFFIX
exit
fi

for ((i=1;i<=$1;i+=1)); do
#echo chunk_$i.$2 
if [ ! -f chunk_$i.$2 ]; then
echo chunk_$i.$2 
fi
done
