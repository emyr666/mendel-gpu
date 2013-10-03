#!/bin/bash

if [ $# -lt 1 ] ; then
  echo "<study>"
  exit 1
fi

#chrs='1 2 3 4 5'
chrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X'
study=$1
chunks=$2

for chr in $chrs
do
  merge_phased.sh $study $chr 1>merge.$chr.out 2>merge.$chr.err &  
done
