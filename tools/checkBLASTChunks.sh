#!/bin/bash
if [ $# -eq 0 ]; then
echo usage: checkChunks NUM SUFFIX EXPECT
exit
fi

for ((i=1;i<=$1;i+=1)); do
if [ -f chunk_$i.$2 ]; then
#echo chunk_$i.$2 
count=`zcat -f chunk_$i.$2  | grep -c "^# Query:"`
if [ $count -ne $3 ]; then
echo chunk_$i.$2
echo $count
fi
fi
done
