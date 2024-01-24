#!/bin/bash
set -e -o pipefail
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
main() {

  echo $(locale)

  #################
  # Download data #
  #################

  echo ""
  echo "=== Setup ==="
  echo "  [*] Downloading input files..."
  dx-download-all-inputs --parallel > /dev/null

  ################
  # Housekeeping #
  ################

  echo "  [*] Performing some housekeeping..."

  # Path setup

  export PATH=$PATH:/stjude/bin
  export PERL5LIB=/stjude/lib/perl:/stjude/bin

  # Cloud config

  export SJ_CONFIGS=/stjude/configs/

  # Variables

  export LEN=$(getReadLength.sh $star_mapped_sorted_bam_path)
  export SHELL=/bin/bash
  export OUT_DIR=/home/dnanexus/out && mkdir -p $OUT_DIR

  #. import_config.sh genome $ref_name 

  ##############################
  # Generate any missing files #
  ##############################

  echo "  [*] Downloading reference files ..."
  local_reference_dir=/stjude/reference
  if [ ! -d $local_reference_dir ]
  then 
    mkdir $local_reference_dir
  fi
  species="Homo_sapiens"
  #ref_name="GRCh37-lite"
  echo "      - Species: $species"
  echo "      - Genome build: $ref_name"
  ref_dir=$local_reference_dir/$species/$ref_name
  if [ ! -d $ref_dir ] 
  then 
    mkdir -p $ref_dir
  fi

  # Download global reference data
  dx download -o $ref_dir -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/$species/$ref_name/2BIT
  dx download -o $ref_dir -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/$species/$ref_name/CICERO
  mkdir -p $ref_dir/CLINICAL/MEDAL_CONFIG
  dx download -o $ref_dir/CLINICAL -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/$species/$ref_name/CLINICAL/GoldGene.lst
  cp $ref_dir/CLINICAL/GoldGene.lst $ref_dir/CLINICAL/MEDAL_CONFIG/
  dx download -o $ref_dir -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/$species/$ref_name/FASTA
  mkdir -p $ref_dir/SUPPORT
  dx download -o $ref_dir/SUPPORT -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/$species/$ref_name/SUPPORT/chromosome_lengths.txt

  # Download CICERO specific reference data
  dx download -o $ref_dir -r project-F5444K89PZxXjBqVJ3Pp79B4:/pipeline/CICERO/reference/$species/$ref_name/mRNA

  # Check for bam index
  if [ ! -d "/home/dnanexus/in/bam_index" ]; then
    echo "  [*] Indexing BAM file..."
    sleep 1 # for pretty output
    export bam_index_path="/home/dnanexus/in/bam_index/$star_mapped_sorted_bam_prefix.bam.bai"
    export bam_index_prefix="$star_mapped_sorted_bam_prefix.bam.bai"
    mkdir -p /home/dnanexus/in/bam_index/
    ( samtools index $star_mapped_sorted_bam_path $bam_index_path )&
    export bam_index_pid=$!
  fi

  if [ -n "$bam_index_pid" ]; then
    wait $bam_index_pid
  fi

  # Linking files where cloud configs expects them

  ln -s $bam_index_path $star_mapped_sorted_bam_path.bai

  if [ -f "$rnapeg_junctions_path" ]; then
    POSSIBLE_JUNCTION_OPTION=" -j $rnapeg_junctions_path "
    echo "  [*] Junctions file: $POSSIBLE_JUNCTION_OPTION"
  else
    POSSIBLE_JUNCTION_OPTION=""
    echo "  [*] Junctions file not supplied"
  fi

  # Debug

  echo "  [*] Summary"
  echo "      - Read LEN: $LEN"
  echo "      - Output directory: $OUT_DIR"
  echo ""

  echo "=== Cicero ==="
  echo "      - Optimize: $optimize"
  echo "      - Soft clip cutoff: $sc_cutoff"
  echo "      - Soft clip shift: $sc_shift"
  echo "      - Reference: $ref_name"
  echo "      - Disable excludes list: $disable_excludes"

  optional_args=
  if [ "$optimize" = "true" ]
  then
    optional_args="-p"
  else
    optional_args="-no-optimize"
    if [ $sc_shift ]
    then
      optional_args="$optional_args -c $sc_shift"
    elif [ $sc_cutoff ]
    then
      optional_args="$optional_args -s $sc__cutoff"
    fi
  fi 
  if [ "$disable_excludes" = "true" ]
  then
    optional_args="$optional_args -d"
  fi
  echo "  [*] Arguments ..."
  echo "      - -n 16 -b $star_mapped_sorted_bam_path -g $ref_name -r /opt/cicero/reference $POSSIBLE_JUNCTION_OPTION -o $OUT_DIR $optional_args"

  docker run -v /home/dnanexus/in:/home/dnanexus/in -v /home/dnanexus/out:/home/dnanexus/out -v /stjude/reference:/opt/cicero/reference ghcr.io/stjude/cicero:CICERO_VERSION Cicero.sh -n 16 -b $star_mapped_sorted_bam_path -g $ref_name -r /opt/cicero/reference $POSSIBLE_JUNCTION_OPTION -o $OUT_DIR $optional_args

  echo ""
  echo "=== Output ==="

  ##############
  # Zip output #
  ##############

  echo "  [*] Zipping up..."
  tar -zcvf cicero_output.tar.gz out/

  cp $OUT_DIR/CICERO_DATADIR/*/final* . 

  ################
  # Upload files #
  ################

  # The following line(s) use the utility dx-jobutil-add-output to format and
  # add output variables to your job's output as appropriate for the output
  # class.  Run "dx-jobutil-add-output -h" for more information on what it
  # does.

  echo "  [*] Uploading files..."
  cicero_output=$(dx upload cicero_output.tar.gz --brief --path "$star_mapped_sorted_bam_prefix".cicero_output.tar.gz)
  cicero_final_fusions=$(dx upload final_fusions.txt --brief --path "$star_mapped_sorted_bam_prefix".final_fusions.txt)
  cicero_final_fusions_html=$(dx upload final_fusions.report.html --brief --path "$star_mapped_sorted_bam_prefix".final_fusions.html)

  dx-jobutil-add-output cicero_output "$cicero_output" --class=file
  dx-jobutil-add-output cicero_final_fusions "$cicero_final_fusions" --class=file
  dx-jobutil-add-output cicero_final_fusions_html "$cicero_final_fusions_html" --class=file

}
# vim: set et ts=2 sts=2 sw=2:
