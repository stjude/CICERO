#!/bin/bash
# QCs the soft clips results for a run
#
# $1 = case EBN 
# $2 = data directory
source qclib.sh

CASE=$1
DIR=$2 

CaseDir="$DIR/$CASE"

starttestcase ITDSClips Using dir $DIR
# Do tests

starttest DirExists
if [ -d "$CaseDir" ]; then passtest; else aborttestcase; fi

cd $CaseDir

starttest AnySclipFilesExist
if ls $CaseDir/*local.SC >&2; then passtest; else aborttestcase; fi

starttest CoverSize
cutoff=200
while read file
do 
  if [ "`filesize $file`" -lt $cutoff ]
  then 
    ls -l $file >&2
    failtestifactive Found at least one cover that was too small
  fi
done< <(ls $CaseDir/*local.SC | sed '$d')
passtestbydefault

summarize
