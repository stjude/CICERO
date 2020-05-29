#!/bin/bash
# Performs setup tasks for running Cicero. 
#
# Accepts the following parameters: 
# $1 = type
# $2 = genome
# $3 = analysis configuration file
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
if [ `which_config.sh app cicero` ]
then
  . import_config.sh app cicero
fi 
if [ `which_config.sh target $TARGET` ]
then 
  . import_config.sh target $TARGET
fi
. import_config.sh genome $GENOME

. steplib.sh

set_step_script_dir $RUN_DIR 

#
# Step 01: extractSClips 
#
init_step extractSClips 
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
touch `get_step_cmds_file`
cat /dev/null > `get_step_cmds_file`

while read case_bam
 do
   bam="$DATA_DIR/\$case_bam/\$case_bam.bam"
   LEN=\`getReadLength.sh \$bam\` 
   get_sc_cmds.pl -i \$bam -o $DATA_DIR/\$case_bam -genome $GENOME -l \$LEN >> `get_step_cmds_file`
   echo "get_geneInfo.pl -i \$bam -o $DATA_DIR/\$case_bam -l \$LEN -genome $GENOME -s \$case_bam " >> `get_step_cmds_file`
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
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$case_bam qc_cicero_sclips.sh \$case_bam  $DATA_DIR 
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

#TODO: Add a step that concatenates the two output sets together
#get_geneInfo.pl and extract_range.pl - DONE
cat > `get_step_local_work_script` <<EOF
#!/bin/bash
while read case_bam 
 do
   bam="$DATA_DIR/\$case_bam/\$case_bam.bam"
   LEN=\`getReadLength.sh \$bam\` 

   prepareCiceroInput.pl -o $DATA_DIR/\$case_bam -genome $GENOME -s 250 -p \$case_bam -l \$LEN -f $DATA_DIR/\$case_bam/\$case_bam.gene_info.txt
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
   
   get_cicero_cmds.pl -i \$bam -genome $GENOME -l \$LEN -o $DATA_DIR/\$case_bam >> `get_step_cmds_file`
#   echo "Cicero.pl -d \$bam -o $DATA_DIR/\$case_bam -l \$LEN" --ref_genome $GENOME -f $DATA_DIR/\$case_bam/\$case_bam.sclip.txt >> `get_step_cmds_file`
 done < $RUN_DIR/config.txt 
EOF
write_step_submit_script

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
# Step 03: annotate
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


cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
touch `get_step_cmds_file`
cat /dev/null > `get_step_cmds_file`

while read case_bam 
 do
   bam="$DATA_DIR/\$case_bam/\$case_bam.bam"
   LEN=\`getReadLength.sh \$bam\` 
   ln -s `which blacklist.genes.txt` $DATA_DIR/\$case_bam
   echo "annotate.pl -i \$bam -o $DATA_DIR/\$case_bam -l \$LEN -genome $GENOME -s \$case_bam -f $DATA_DIR/\$case_bam/\$case_bam.gene_info.txt -j $DATA_DIR/\$case_bam/\$case_bam.bam.junctions.tab.shifted.tab" >> `get_step_cmds_file`
   echo "annotate.pl -i \$bam -o $DATA_DIR/\$case_bam -l \$LEN -genome $GENOME -s \$case_bam -f $DATA_DIR/\$case_bam/\$case_bam.gene_info.txt -internal" >> `get_step_cmds_file`
 done < $RUN_DIR/config.txt 
EOF
write_step_submit_script



#
# Step 04: filter
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

# Step 05
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

cat > `get_step_local_work_script` <<EOF
#!/bin/bash
while read case_bam 
 do
   cat $HTML_FIRST_HALF $DATA_DIR/\$case_bam/final_fusions.txt $HTML_SECOND_HALF> $DATA_DIR/\$case_bam/final_fusions.report.html 
 done < $RUN_DIR/config.txt 
EOF


