#!/bin/bash
# QCs the annotation results for a run
#
# $1 = bam EBN 
# $2 = data directory
source qclib.sh

EBN=$1
DIR=$2 

starttestcase Annotations Using dir $DIR/$EBN
# Do tests

starttest DirExists
if [ -d "$DIR/$EBN" ]; then passtest; else aborttestcase; fi

cd $DIR/$EBN

starttest AnyFilesExist
if ls $DIR/$EBN/annotated.*.txt >&2; then passtest; else aborttestcase; fi

starttest AnnotationSize 
cutoff=0
while read file
do 
  if [ "`filesize $file`" -lt $cutoff ]
  then 
    ls -l $file >&2
    failtestifactive Found at least one cover that was too small
  fi
done< <(ls $DIR/$EBN/annotated.*.txt | sed '$d' )
passtestbydefault

summarize
