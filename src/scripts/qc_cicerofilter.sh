#!/bin/bash
# QCs the annotation results for a run
#
# $1 = data directory
# $2 = bam EBN
source qclib.sh

DIR=$1
EBN=$2

starttestcase Filter Using dir $DIR
# Do tests

starttest DirExists
if [ -d "$DIR" ]; then passtest; else aborttestcase; fi

cd $DIR

starttest AnyFilesExist
if ls $DIR/$EBN/final* >&2; then passtest; else aborttestcase; fi

starttest ResultSize
cutoff=200
while read file
do
  if [ "`filesize $file`" -lt $cutoff -a ! -L $file ]
  then
    ls -l $file >&2
    failtestifactive Results file $file too small
  fi
done< <(find $DIR/$EBN -type f -name "final*")
passtestbydefault

summarize
