#!/bin/bash
# Performs setup tasks for running Cicero ITD.
#
# Directory structure:
# <data directory>
# +-- <sample>
#     +-- <sample>.bam (input)
#     +-- <sample>.bam.bai (input)
#     +-- (outputs written here)
#
# The run config file is simply a list of samples, one per line
#
# Accepts the following parameters: 
# $1 = type
# $2 = genome
# $3 = run configuration file
# $4 = run directory
# $5 = data directory

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
ARGS=
while [ "$#" -gt 0 ]
do
  case "$1" in
  -*)             echo "Unrecognized option: $1" >&2; exit 1 ;;
  *)              ARGS="$ARGS $1"
  esac
  shift
done
set -- $ARGS
TARGET=$1
GENOME=$2
ANLS_CONFIG=$3
RUN_DIR=$4
DATA_DIR=$5

cp $ANLS_CONFIG $RUN_DIR/config.txt

# Source the 3 relevant config files in order.
echo
echo "Reading config files (you may see some error messages--these are OK):"
echo "* Application level"
echo "* Sequencing type level"
echo "..."
if [ `which_config.sh app cicero-itd` ]
then
  . import_config.sh app cicero-itd
fi 
if [ `which_config.sh target $TARGET` ]
then 
  . import_config.sh target $TARGET
fi
. import_config.sh genome $GENOME

# Cicero ITD uses a custom gene list and gene model.  These will be generated in
# the first step and reused later.
GENE_MODEL_DIR=$RUN_DIR/gene_model

. steplib.sh

set_step_script_dir $RUN_DIR 

init_step geneInfo

cat > `get_step_local_work_script` <<EOF
#!/bin/bash

# Make the temporary gene model
# NOTE: this used to be KNOWN_ITD_FILE instead of CLINCLS_GOLD_GENE_LIST_FILE
mkdir -p $GENE_MODEL_DIR
cp $CLINCLS_GOLD_GENE_LIST_FILE $GENE_MODEL_DIR/genes.lst
cat $REFSEQ_REFFLAT | grep -w -f $CLINCLS_GOLD_GENE_LIST_FILE > $GENE_MODEL_DIR/tmp_gene_model

exit 0
EOF

cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
touch `get_step_cmds_file`
cat /dev/null > `get_step_cmds_file`

while read case_bam
 do
   bam="$DATA_DIR/\$case_bam/\$case_bam.bam"
   LEN=\`getReadLength.sh \$bam\` 
   echo "get_geneInfo.pl -i \$bam -o $DATA_DIR/\$case_bam -l \$LEN -genome $GENOME -p \$case_bam -genemodel $GENE_MODEL_DIR/tmp_gene_model" >> `get_step_cmds_file`
 done < $RUN_DIR/config.txt
EOF
write_step_submit_script

#
# Step 01: extractSClips 
#
init_step extractSClips 
cat > `get_step_qc_script` <<EOF
#!/bin/bash
# QC:
anyfail=no
while read case_bam
do
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$case_bam qc_cicero_gInfo.sh \$case_bam  $DATA_DIR 
  then anyfail=yes
  fi
done < $RUN_DIR/config.txt
if [ "\$anyfail" == "yes" ]
then 
  echo There were QA failures
  echo "Exiting..."
  exit 1
fi
EOF


cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
touch `get_step_cmds_file`
cat /dev/null > `get_step_cmds_file`

while read case_bam
 do
   bam="$DATA_DIR/\$case_bam/\$case_bam.bam"
   LEN=\`getReadLength.sh \$bam\` 
   if [ "$ANALYTE" == "DNA" ]
   then 
     echo "extractSClip.gene.pl -i \$bam -f $DATA_DIR/\$case_bam/\$case_bam.gene_info.txt -o $DATA_DIR/\$case_bam -genome $GENOME -l \$LEN -m 2 -DNA" >> `get_step_cmds_file`
   else
     echo "extractSClip.gene.pl -i \$bam -f $DATA_DIR/\$case_bam/\$case_bam.gene_info.txt -o $DATA_DIR/\$case_bam -genome $GENOME -l \$LEN -m 2" >> `get_step_cmds_file`
   fi 
 done < $RUN_DIR/config.txt
EOF
write_step_submit_script

#
# Step 02: run Cicero
#
init_step cicero

cat > `get_step_qc_script` <<EOF
#!/bin/bash
# QC:
anyfail=no
while read case_bam
do
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$case_bam qc_cicero_itd_sclips.sh \$case_bam  $DATA_DIR 
  then anyfail=yes
  fi
done < $RUN_DIR/config.txt
if [ "\$anyfail" == "yes" ]
then 
  echo There were QA failures
  echo "Exiting..."
  exit 1
fi
EOF

cat > `get_step_local_work_script` <<EOF
#!/bin/bash

THRESHOLD=200000
SC_CUTOFF=3

while read case_bam
  do
    bam="$DATA_DIR/\$case_bam/\$case_bam.bam"
    LEN=\`getReadLength.sh \$bam\`
    SOFTCLIP_COUNT=\`wc -l $DATA_DIR/\$case_bam/*.cover | tail -n 1 | awk '{print \$1}'\`

    sc_cutoff_arg=
    if [ \$THRESHOLD -gt 0 ]
    then
       if [ \$SOFTCLIP_COUNT -gt \$THRESHOLD ]
       then
          sc_cutoff_arg="-m \$SC_CUTOFF"
       fi
    fi

    prepareCiceroInput.pl -o $DATA_DIR/\$case_bam -genome $GENOME -s 250 -p local -l \$LEN -f $DATA_DIR/\$case_bam/\$case_bam.gene_info.txt \$sc_cutoff_arg
  done < $RUN_DIR/config.txt
EOF

cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
touch `get_step_cmds_file`
cat /dev/null > `get_step_cmds_file`

while read case_bam 
 do
   bam="$DATA_DIR/\$case_bam/\$case_bam.bam"
   LEN=\`getReadLength.sh \$bam\` 
   
   get_cicero_cmds.pl -i \$bam -genome $GENOME -l \$LEN -o $DATA_DIR/\$case_bam -c 10 -out_prefix=local >> `get_step_cmds_file`
 done < $RUN_DIR/config.txt 
EOF
write_step_submit_script

#
# Step 03: combine
#

init_step combine
cat > `get_step_local_work_script` <<EOF
#!/bin/bash
while read case_bam
  do
    bam="$DATA_DIR/\$case_bam/\$case_bam.bam"
    cat $DATA_DIR/\$case_bam/*/unfiltered.fusion.txt > $DATA_DIR/\$case_bam/unfiltered.fusion.txt
    cat $DATA_DIR/\$case_bam/*/unfiltered.internal.txt > $DATA_DIR/\$case_bam/unfiltered.internal.txt
  done < $RUN_DIR/config.txt
EOF

#
# Step 04: annotate
#
init_step annotate

cat > `get_step_qc_script` <<EOF
#!/bin/bash
# QC:
anyfail=no
while read case_bam 
do
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$case_bam qc_cicero.sh \$case_bam $DATA_DIR 
  then anyfail=yes
  fi
done < $RUN_DIR/config.txt
if [ "\$anyfail" == "yes" ]
then 
  echo There were QA failures
  echo "Exiting..."
  exit 1
fi
EOF

echo "ANALYTE is $ANALYTE"
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
touch `get_step_cmds_file`
cat /dev/null > `get_step_cmds_file`

while read case_bam 
 do
   bam="$DATA_DIR/\$case_bam/\$case_bam.bam"
   LEN=\`getReadLength.sh \$bam\` 
   ln -s $EXCLUDED_GENES $DATA_DIR/\$case_bam
   if [ "$ANALYTE" == "DNA" ]
   then 
     echo "annotate.pl -i \$bam -o $DATA_DIR/\$case_bam -f $DATA_DIR/\$case_bam/\$case_bam.gene_info.txt -l \$LEN -genome $GENOME -s \$case_bam -internal -DNA" >> `get_step_cmds_file`
     echo "annotate.pl -i \$bam -o $DATA_DIR/\$case_bam -f $DATA_DIR/\$case_bam/\$case_bam.gene_info.txt -l \$LEN -genome $GENOME -s \$case_bam" >> `get_step_cmds_file`
   else
     echo "annotate.pl -i \$bam -o $DATA_DIR/\$case_bam -f $DATA_DIR/\$case_bam/\$case_bam.gene_info.txt -l \$LEN -genome $GENOME -s \$case_bam" >> `get_step_cmds_file`
     echo "annotate.pl -i \$bam -o $DATA_DIR/\$case_bam -f $DATA_DIR/\$case_bam/\$case_bam.gene_info.txt -l \$LEN -genome $GENOME -s \$case_bam -internal" >> `get_step_cmds_file`
   fi 
 done < $RUN_DIR/config.txt 
EOF
write_step_submit_script



#
# Step 05: filter
#
init_step filter

cat > `get_step_qc_script` <<EOF
#!/bin/bash
# QC:
anyfail=no
while read case_bam 
do
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$case_bam qc_annotation.sh \$case_bam $DATA_DIR 
  then anyfail=yes
  fi
done < $RUN_DIR/config.txt
if [ "\$anyfail" == "yes" ]
then 
  echo There were QA failures
  echo "Exiting..."
  exit 1
fi
EOF


cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
while read case_bam 
do
  echo cicero_filter.sh $DATA_DIR \$case_bam $GENOME
done < $RUN_DIR/config.txt > `get_step_cmds_file`
EOF
write_step_submit_script

#
# Step 06
#
init_step final_qa

cat > `get_step_qc_script` <<EOF
#!/bin/bash
# QC:
echo -n "" > final_qa.txt
anyfail=no
while read case_bam
do
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$case_bam qc_cicerofilter.sh $DATA_DIR/ \$case_bam
  then 
    anyfail=yes
    echo "FAIL \$case_bam" >> final_qa.txt
  else
    echo "PASS \$case_bam" >> final_qa.txt 
  fi
done < $RUN_DIR/config.txt
if [ "\$anyfail" == "yes" ]
then 
  echo There were QA failures
  echo "Exiting..."
  exit 1
fi
EOF

