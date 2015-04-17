#!/bin/bash
if [ -f $1 ]
then
    source $1
    exit 0
fi

for i in {1..10}
do
  if [[ "$i" < 6 ]]
  then
    sleep 30
  else
    sleep 180
  fi
  if [ -f $1 ]
  then
    source $1
    exit 0
  fi
done

exit 1

