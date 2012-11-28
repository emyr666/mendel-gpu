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
  chroms='X'
  #chroms='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X'
elif [ $useX -eq 0 ] ; then
  chroms='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22'
else
  echo use X is 1 or 0
  exit 1
fi
chunksize=1000
regionsize=10000

rm -f $IMPUTATION_DATA'/'$dbname'/permutation'*
rm -f $IMPUTATION_DATA'/'$dbname'/sex'*
mysql -h $DBHOST -u $DBUSER --password=$DBPW $dbname << END > $IMPUTATION_DATA'/'$dbname'/permutation'
select seq,perm from person order by seq;
END
mysql -h $DBHOST -u $DBUSER --password=$DBPW $dbname << END > $IMPUTATION_DATA'/'$dbname'/sex'
select seq,sex from person order by perm;
END
#exit 0
for chr in $chroms
do
  fetchdata.sh $dbname $chr $regionsize $chunksize
done
