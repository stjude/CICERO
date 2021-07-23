#!/bin/bash
# QCs the soft clips results for a run
#
# $1 = case EBN
# $2 = data directory
source qclib.sh

CASE=$1
DIR=$2

CaseDir="$DIR/$CASE"

starttestcase SClips Using dir $DIR
# Do tests

starttest DirExists
if [ -d "$CaseDir" ]; then passtest; else aborttestcase; fi

cd $CaseDir

starttest AnySclipFilesExist
if ls $CaseDir/$CASE.*.cover >&2; then passtest; else aborttestcase; fi

starttest AnyGinfoFilesExist
if ls $CaseDir/$CASE.gene_info.txt >&2; then passtest; else aborttestcase; fi

starttest GinfoSize
cutoff=200
while read file
do
  if [ "`filesize $file`" -lt $cutoff ]
  then
    ls -l $file >&2
    failtestifactive Found at least one cover that was too small
  fi
done< <(ls $CaseDir/$CASE.gene_info.txt | sed '$d' )
passtestbydefault


starttest CoverSize
anynonempty=
# anysmall=
cutoff=200
while read file
do
  if [ -s $file ]
  then anynonempty=1
  fi
  ## B2P - Removed as it was causing false problems
  ## (Note that we already know that the coverage file is non-empty,
  ## and there is no header information)
  # if [ "`filesize $file`" -lt $cutoff ]
  # then
  #   anysmall=1
  #   ls -l $file >&2
  # fi
done< <(ls $CaseDir/$CASE.*.cover | sed '$d')
if [ ! "$anynonempty" ]
then failtest All cover files were empty
# elif [ "$anysmall" ]
# then warntest Found at least one cover that was too small
else passtest
fi

summarize
