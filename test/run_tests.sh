#!/usr/bin/env bash

##########################
### Verify environment ###
##########################
export PATH=`readlink -f $(pwd)/../src/scripts`:$PATH
export PERL5LIB=`readlink -f $(pwd)/../src/perllib`:`readlink -f $(pwd)/../dependencies/lib/perl`:$PERL5LIB
export CLASSPATH=`readlink -f $(pwd)/../src/javalib`/*
export LC_ALL=C

##################
### Set inputs ###
##################
SAMPLE=SJBALL020016_C27.1.1
BAMFILE=`readlink -f $(pwd)/data/input/$SAMPLE/${SAMPLE}.bam`
REFDIR=`readlink -f $(pwd)/../reference`
GENOME=GRCh37-lite
JUNCTIONS=`readlink -f $(pwd)/data/input/$SAMPLE/${SAMPLE}.bam.junctions.tab`

####################
### Run the code ###
####################
rm -rf ${SAMPLE}.tmp
time Cicero.sh -b $BAMFILE -r $REFDIR -g $GENOME -n 4 -j $JUNCTIONS -o ${SAMPLE}.$$.tmp

########################
### Test the results ###
########################
tmp1=$(mktemp)
tmp2=$(mktemp)

RETURN_CODE=0
for GS_FILENAME in data/output/$SAMPLE/*; do
    if [[ ! -f $GS_FILENAME ]]; then continue; fi
    FILENAME=$(basename $GS_FILENAME)
    AT_FILENAME=${SAMPLE}.$$.tmp/CICERO_DATADIR/$SAMPLE/$FILENAME

    if [[ ! -f $AT_FILENAME ]]; then
        echo "DNE .... $FILENAME"
        RETURN_CODE=1
        continue
    fi

    sort $GS_FILENAME > $tmp1
    sort $AT_FILENAME > $tmp2
    cmp $tmp1 $tmp2 > /dev/null 2>&1
    if [[ $? == 0 ]]; then
        echo "ok ..... $FILENAME"
    else
        echo "FAIL ... $FILENAME"
        echo $GS_FILENAME $AT_FILENAME
        RETURN_CODE=1
    fi
done

rm $tmp1
rm $tmp2

exit $RETURN_CODE
