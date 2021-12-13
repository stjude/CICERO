#!/bin/bash
# QCs the soft clips results for a run
#
# $1 = case EBN
# $2 = data directory
source qclib.sh

CASE=$1
DIR=$2

CaseDir="$DIR/$CASE"

starttestcase gInfo Using dir $DIR
# Do tests

starttest DirExists
if [ -d "$CaseDir" ]; then passtest; else aborttestcase; fi

cd $CaseDir

starttest AnyGinfoFilesExist
if ls $CaseDir/$CASE.gene_info.txt >&2; then passtest; else aborttestcase; fi

starttest GinfoSize
cutoff=200
while read file
do
  if [ "`filesize $file`" -lt $cutoff ]
  then
    ls -l $file >&2
    failtestifactive Found at least one gene_info.txt file that was too small
  fi
done< <(ls $CaseDir/$CASE.gene_info.txt | sed '$d' )
passtestbydefault

summarize
