#!/bin/bash
# Encapsulates the multiple steps that filtering and annotating Cicero 
# output has become
# 
# Parameters: 
# DATA_DIR  Location of the data files from cicero 
# CASE_BAM  Name of the sample to run on 
# GENOME    Genome of the sample being run

# Show usage information if no parameters were sent
if [ "$#" -lt 3 ]; then about.sh $0; exit 1; fi

DATA_DIR=$1
CASE_BAM=$2
GENOME=$3

cat $DATA_DIR/$CASE_BAM/${CASE_BAM}_annotated.fusion.txt $DATA_DIR/$CASE_BAM/${CASE_BAM}_annotated.internal.txt > $DATA_DIR/$CASE_BAM/${CASE_BAM}_annotated.all.txt
sv_inframe.pl -genome $GENOME -single $DATA_DIR/$CASE_BAM/${CASE_BAM}_annotated.all.txt -fq
bam="$DATA_DIR/$CASE_BAM/$CASE_BAM.bam"
LEN=`getReadLength.sh $bam`
rank_SVs.pl -i $DATA_DIR/$CASE_BAM/${CASE_BAM}_annotated.all.txt.frame.tab -genome $GENOME -l $LEN

# rename to get the sample name into the file name.
mv $DATA_DIR/$CASE_BAM/final_fusions.txt $DATA_DIR/$CASE_BAM/${CASE_BAM}_final_fusions.txt

head -n 1 $DATA_DIR/$CASE_BAM/${CASE_BAM}_final_fusions.txt > $DATA_DIR/$CASE_BAM/${CASE_BAM}_final_internal.txt
