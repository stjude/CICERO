#!/bin/bash
# Writes scripts for Cicero post-proces runs.
#
# Usage: setup_cicero_post.sh [OPTIONS] TARGET PROJECT CONFIG RUN_DIR DATA_DIR
#
# Parameters:
#
# $1 = target
# $2 = genome
# $3 = analysis configuration file
# $4 = run directory
# $5 = data directory

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
TARGET=$1
GENOME=$2
ANLS_CONFIG=$3
RUN_DIR=$4
DATA_DIR=$5

# Source the 3 relevant config files in order.
echo
echo "Reading config files (you may see some error messages--these are OK):"
echo "* Application level"
echo "* Sequencing target level"
echo "..."
. import_config.sh app cicero-post 
. import_config.sh target $TARGET
. import_config.sh genome $GENOME

# Read step script
. steplib.sh

# Validate analysis config
if [ ! -f $ANLS_CONFIG ]
then echo "Analysis config not found: $ANLS_CONFIG" >&2; exit 1
fi

# Get config name and containing dir
ANLS_CONFIG_PATH=`readlink -f $ANLS_CONFIG`
ANLS_CONFIG_NAME=`basename $ANLS_CONFIG`
ANLS_CONFIG_DIR=`dirname $ANLS_CONFIG_PATH`
ANLS_CONFIG_DIR_NAME=`basename $ANLS_CONFIG_DIR`

# Set the run directory
set_step_script_dir $RUN_DIR

init_step fusion_builder
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
# This script writes fusion builder commands to a file. 

# Clean up the commands file
cat /dev/null > `get_step_cmds_file`

while read case_bam
do 
  bam="$DATA_DIR/\$case_bam/\$case_bam.bam"

  echo "build_fusions.sh $GENOME $DATA_DIR/\$case_bam/final_fusions- -i $DATA_DIR//\$case_bam/final_fusions.txt -h -x -S \$case_bam " >> `get_step_cmds_file`
  echo "build_fusions.sh $GENOME $DATA_DIR/\$case_bam/final_internal- -i $DATA_DIR//\$case_bam/final_internal.txt -h -x -S \$case_bam " >> `get_step_cmds_file`

done < $ANLS_CONFIG
EOF
write_step_submit_script


init_step append-allele
cat > `get_step_qc_script` <<EOF
#!/bin/bash
# This QC's the output of fusion builder.
anyfail=no
while read case_bam
do
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$case_bam qc_fusions.sh $DATA_DIR/\$case_bam final_fusions final_fusions
  then 
    anyfail=yes
  fi
done < $ANLS_CONFIG
if [ "\$anyfail" == "yes" ]
then 
  echo There were QA failures
  echo "Exiting..."
  exit 1
fi
EOF

cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
# This script writes append allele count commands to a file. 

# Clean up the commands file
cat /dev/null > `get_step_cmds_file`

while read case_bam 
do 
  bam="$DATA_DIR/\$case_bam/\$case_bam.bam"

  echo "java.sh org.stjude.compbio.sv.counting.AppendAlleleCounts -i $DATA_DIR/\$case_bam/final_fusions-event_fusion.txt -b \$bam -o $DATA_DIR/\$case_bam/final_fusions.counts -f $FASTA -V SILENT -l 20" >> `get_step_cmds_file`
  echo "java.sh org.stjude.compbio.sv.counting.AppendAlleleCounts -i $DATA_DIR/\$case_bam/final_internal-event_fusion.txt -b \$bam -o $DATA_DIR/\$case_bam/final_internal.counts -f $FASTA  -V SILENT -l 20" >> `get_step_cmds_file`

done < $ANLS_CONFIG
EOF
write_step_submit_script


init_step classification
cat > `get_step_qc_script` <<EOF
#!/bin/bash
# This QC's the output of append allele count.
anyfail=no
while read case_bam 
do
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$case_bam qc_append_allele.sh $DATA_DIR/\$case_bam final_fusions final_fusions
  then 
    anyfail=yes
  fi
done < $ANLS_CONFIG
if [ "\$anyfail" == "yes" ]
then 
  echo There were QA failures
  echo "Exiting..."
  exit 1
fi
EOF

cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
# This script writes medal ceremony commands to a file. 

# Clean up the commands file
cat /dev/null > `get_step_cmds_file`

while read case_bam
do 
  bam="$DATA_DIR/\$case_bam/\$case_bam.bam"

  echo "cd $DATA_DIR/\$case_bam/; medal_ceremony_using_configs.sh GRCh37-lite -outfile-fq -single-sv $DATA_DIR/\$case_bam/final_fusions.counts" >> `get_step_cmds_file`
  echo "cd $DATA_DIR/\$case_bam/; medal_ceremony_using_configs.sh GRCh37-lite -outfile-fq -single-sv $DATA_DIR/\$case_bam/final_internal.counts" >> `get_step_cmds_file`

done < $ANLS_CONFIG
EOF
write_step_submit_script

init_step final_qa

cat > `get_step_qc_script` <<EOF
#!/bin/bash
# QC:
echo -n "" > $RUN_DIR/final_qa.txt
anyfail=no
while read case_bam 
do
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$case_bam qc_classification.sh $DATA_DIR/\$case_bam final_fusions final_fusions 
  then 
    anyfail=yes
    echo "FAIL \$case_bam" >> $RUN_DIR/final_qa.txt
  else
    echo "PASS \$case_bam" >> $RUN_DIR/final_qa.txt 
  fi
done < $ANLS_CONFIG
if [ "\$anyfail" == "yes" ]
then 
  echo There were QA failures
  echo "Exiting..."
  exit 1
fi
EOF
