#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "Usage <dbname> <use X?[1|0]>"
  exit 1
fi

dbname=$1
useX=$2
mkdir -p $IMPUTATION_DATA'/'$dbname
ln -fs $IMPUTATION_DATA'/'$dbname

if [ $useX -eq 1 ] ; then
  chroms='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X'
elif [ $useX -eq 0 ] ; then
  chroms='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22'
else
  echo use X is 1 or 0
  exit 1
fi

for chr in $chroms
do
  fetchdata.sh $dbname $chr 1>fetch.$chr.out 2>fetch.$chr.err &  
done
