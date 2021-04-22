#!/usr/bin/env bash

exit_trap() {
    for PID in $(pgrep -P $$); do
        PS_LINE=$(ps --no-header -p $PID)
        if [[ $PS_LINE ]]; then
            echo "Killing $PS_LINE"
            kill -9 $PID
        fi
    done
}
trap exit_trap EXIT

MAX_PORT_NUM=32000
MIN_PORT_NUM=2000

BAMFILE=""
JUNCTIONS=""
REFFA=""
GENOME=""
SAMPLE=""
TARGET="TRANSCRIPTOME"
NCORES=
OUTDIR=$(pwd)
THRESHOLD=200000
SC_CUTOFF=3
SC_SHIFT=
OPTIMIZE=1
DISABLE_EXCLUDE=0

usage() {
    echo "Cicero.sh [-h] [-n ncores] -b bamfile -g genome -r refdir [-j junctions] [-o outdir] [-t threshold] [-s sc_cutoff] [-c sc_shift] [-p]"
    echo "-p - optimize CICERO, sets sc_cutoff=3 and sc_shift=10 [default true]" 
    echo "-s <num> - minimum number of soft clip support required [default=2]"
    echo "-t <num> - threshold for enabling increased soft clip cutoff [default=200000]"
    echo "-c <num> - clustering distance for grouping similar sites [default=3]"
    echo "-j <file> - junctions file from RNApeg"
    echo "-n <num> - number of cores to utilize with GNU parallel"
    echo "-d - disable excluded regions file use"
}


###################
### Read inputs ###
###################
while [ ! -z "$1" ]; do
    case "$1" in
        -h) usage && exit 0; shift;;
        -n) NCORES=$2; shift;;
        -b) BAMFILE=$2; shift;;
        -j) JUNCTIONS=$2; shift;;
        -r) REFDIR=$2; shift;;
        -g) GENOME=$2; shift;;
        -o) OUTDIR=$2; shift;;
        -t) THRESHOLD=$2; shift;;
        -s) SC_CUTOFF=$2; shift;;
        -c) SC_SHIFT=$2; shift;;
        -p) OPTIMIZE=1;;
        -no-optimize) OPTIMIZE=0;;
        -d) DISABLE_EXCLUDE=1;;
    esac
    shift
done

PARALLEL_ARG=
if [ $NCORES ]
then
  PARALLEL_ARG="-j $NCORES"
fi

bam_dir=$(dirname ${BAMFILE})
bam_name=$(basename ${BAMFILE} ".bam")
if [[ -f "${bam_dir}/${bam_name}.bai" ]]
then
  ln -s ${bam_dir}/${bam_name}.bai ${BAMFILE}.bai
fi

#######################
### Validate inputs ###
#######################
if [[ ! -f $BAMFILE ]]; then
    >&2 echo "ERROR: Bamfile '$BAMFILE' does not exist"
    usage
    exit 1
elif [[ ! -f $BAMFILE.bai ]]; then
    >&2 echo "ERROR: Bam index (.bai) '$BAMFILE.bai' does not exist"
    usage
    exit 1
elif [[ ! -d $REFDIR ]]; then
    >&2 echo "ERROR: Reference directory '$REFFA' does not exist"
    usage
    exit 1
elif [[ ! $GENOME ]]; then
    >&2 echo "ERROR: GENOME not defined"
    usage
    exit 1
elif [[ $JUNCTIONS && ! -f $JUNCTIONS ]]; then
    >&2 echo "ERROR: Junctions file '$JUNCTIONS' does not exist"
    usage
    exit 1
fi

SAMPLE=$(basename $BAMFILE .bam)
#if [[ $GENOME != "GRCh37-lite" ]]; then
#    >&2 echo "ERROR: Only genome 'GRCh37-lite' is currently supported"
#    exit 1
#fi

# Check for GNU parallel and BLAT 
which parallel 2> /dev/null > /dev/null
RET=$?
if [ $RET -eq 1 ] 
then
    >&2 echo "ERROR: GNU parallel required"
    exit 1
fi
which gfServer 2> /dev/null > /dev/null
RET=$?
if [ $RET -eq 1 ] 
then
    >&2 echo "ERROR: BLAT required (gfServer)"
    exit 1
fi

cluster_arg=
if [ $SC_SHIFT ] 
then
    echo "Setting SC_SHIFT=$SC_SHIFT"
    cluster_arg="-c $SC_SHIFT"
fi
if [ $OPTIMIZE -eq 1 ]
then
    echo "Optimize: setting SC_SHIFT=10, SC_CUTOFF=3, THRESHOLD=200000"
    SC_SHIFT=10
    SC_CUTOFF=3
    THRESHOLD=200000
    cluster_arg="-c $SC_SHIFT"
fi


# Set inputs
mkdir -p $OUTDIR
cd $OUTDIR
CICERO_RUNDIR=CICERO_RUNDIR
CICERO_DATADIR=CICERO_DATADIR
CICERO_CONFIG=CICERO_CONFIG

##############################
### Construct config files ###
##############################

# Create the config location and set variables to it
CICERO_ROOT=`readlink -f $(dirname ${BASH_SOURCE[0]})/../../`
# TODO - need to totall revamp the config files
export SJ_CONFIGS=configs.${HOSTNAME}.$$.tmp
echo "$CICERO_ROOT/configs"
echo "$SJ_CONFIGS"
cp -rf $CICERO_ROOT/configs $SJ_CONFIGS
echo "SJ_CONFIGS=$SJ_CONFIGS"

# Construct the actual config file
for CONFIG_TYPE_DIR in $SJ_CONFIGS/*; do
    CONFIG_TYPE=$(basename $CONFIG_TYPE_DIR)
    for CONFIG_TEMPLATE in $CONFIG_TYPE_DIR/*.template.txt; do
        CONFIG_FILE=$CONFIG_TYPE_DIR/$(basename $CONFIG_TEMPLATE .template.txt).config.txt
        while read KEY VALUE; do
            if [[ $(echo $VALUE | grep '/') ]]; then
                echo -e "$KEY\t$REFDIR$VALUE"
            elif [[ $KEY == BLAT_PORT && $VALUE == 0 ]]; then
                PORT_NUM=$$
                # Check if the port is open or not in our declared port range, if not, choose a random number until we find a free port
                while [[ $(netstat -tulpn 2> /dev/null | grep LISTEN | grep ":$PORT_NUM\b") || ( $PORT_NUM -le $MIN_PORT_NUM || $PORT_NUM -ge $MAX_PORT_NUM )]]; do
                    echo "WARN: Port $PORT_NUM is not available for use. Trying new port..."
                    # Due to the range of RANDOM, we need to narrow it down using some math
                    RANGE=$(($MAX_PORT_NUM-$MIN_PORT_NUM+1))    # Find the total possible numbers we want
                    RAND_PORT_NUM=$RANDOM                       # Assign out the $RANDOM var 
                    let "RAND_PORT_NUM %= $RANGE"               # Mod the $RANDOM var using the range
                    PORT_NUM=$(($RAND_PORT_NUM+$MIN_PORT_NUM)) # Add the resulting number to the min to get it within the range
                done
                echo -e "$KEY\t$PORT_NUM"
            else
                echo -e "$KEY\t$VALUE"
            fi
        done < $CONFIG_TEMPLATE > $CONFIG_FILE
    done
done

FASTA=$(awk -F$'\t' '$1=="FASTA"{print $2}' $SJ_CONFIGS/genome/${GENOME}.config.txt)
TWOBIT=$(awk -F$'\t' '$1=="TWOBIT"{print $2}' $SJ_CONFIGS/genome/${GENOME}.config.txt)
GFSERVER_LOG=gfServer.$(hostname).$$

BLAT_HOST=$(awk -F$'\t' '$1=="BLAT_HOST"{print $2}' $SJ_CONFIGS/genome/${GENOME}.config.txt)
BLAT_PORT=$(awk -F$'\t' '$1=="BLAT_PORT"{print $2}' $SJ_CONFIGS/genome/${GENOME}.config.txt)

##################
### Start BLAT ###
##################

# Start blat server
gfServer status $BLAT_HOST $BLAT_PORT 1> ${GFSERVER_LOG}.out 2> ${GFSERVER_LOG}.err
RETURN_CODE=$?
if [[ $RETURN_CODE != 0 ]]; then
    gfServer start $BLAT_HOST $BLAT_PORT "$TWOBIT" -log=${GFSERVER_LOG}.log -stepSize=5 1>> ${GFSERVER_LOG}.out 2>> ${GFSERVER_LOG}.err  &
    BLAT_SERVER_PID=$!
    RETURN_CODE=1
    echo "Starting local blat server"
    NUM_CHECKS=0
    while [[ $RETURN_CODE != 0 && $NUM_CHECKS -le 200 ]]; do
        echo "Blat server not yet running ..."
        if  ps -p $BLAT_SERVER_PID > /dev/null; then
            gfServer status $BLAT_HOST $BLAT_PORT 1>> ${GFSERVER_LOG}.out 2>> ${GFSERVER_LOG}.err
            RETURN_CODE=$?
            NUM_CHECKS=$(($NUM_CHECKS + 1))
            sleep 10
        else
            >&2 echo "ERROR: The server has probably died. Please review error logs at below location:"
            >&2 echo "$OUTDIR/${GFSERVER_LOG}.err"
            exit 1
        fi
    done
    if [[ $RETURN_CODE != 0 ]]; then
        >&2 echo "ERROR: Unable to spin up blat server"
        exit 1
    fi
    echo "Blat server is up!"
else
    echo "Blat server is already up!"
fi

# Make config
echo "$SAMPLE" > $CICERO_CONFIG

# Prep datadir
rm -rf $CICERO_DATADIR
mkdir -p $CICERO_DATADIR
mkdir -p $CICERO_DATADIR/$SAMPLE
ln -s $BAMFILE $CICERO_DATADIR/$SAMPLE
ln -s ${BAMFILE}.bai $CICERO_DATADIR/$SAMPLE
ln -s $JUNCTIONS $CICERO_DATADIR/$SAMPLE

# Make rundir
rm -rf $CICERO_RUNDIR
mkdir -p $CICERO_RUNDIR

###############################
### Step 01 - extractSClips ###
###############################

exclude_arg=
if [ $DISABLE_EXCLUDE -eq 1 ]
then
  echo "Disabling excluded regions list"
  exclude_arg="-disable_excludes"
fi

echo "Step 01 - $(date +'%Y.%m.%d %H:%M:%S') - ExtractSClips"
{
LEN=$(getReadLength.sh $BAMFILE)
get_sc_cmds.pl -i $BAMFILE -genome $GENOME -o $CICERO_DATADIR/$SAMPLE -l $LEN $cluster_arg $exclude_arg > cmds-01.sh
echo "get_geneInfo.pl -i $BAMFILE -o $CICERO_DATADIR/$SAMPLE -l $LEN -genome $GENOME -s $SAMPLE" >> cmds-01.sh
parallel --joblog 01_ExtractSClips.log $PARALLEL_ARG < cmds-01.sh
} 1> 01_ExtractSClips.out 2> 01_ExtractSClips.err

## QC
shopt -s nullglob
cover_files=($CICERO_DATADIR/$SAMPLE/*.cover)
cover_file_count=${#cover_files[@]}
softclip_count=$(wc -l  $CICERO_DATADIR/$SAMPLE/*.cover | tail -n 1 | awk '{print $1}')
geneinfo=($CICERO_DATADIR/$SAMPLE/*.gene_info.txt)
geneinfo_count=${#geneinfo[@]}

# Check to ensure the number of output files is equal to the expected number
if [ $(wc -l cmds-01.sh | awk '{ print $1} ') -ne $(( $cover_file_count + $geneinfo_count )) ]
then
  exit "Error in ExtractSClips"
fi

for file in ${cover_files[@]}
do 
   if [ ! -s $file ]
   then
     echo "$file has no soft clipped reads"
   fi
done

sc_cutoff_arg=
if [ $THRESHOLD -gt 0 ]
then
  if [ $softclip_count -gt $THRESHOLD ]
  then
     sc_cutoff_arg="-m $SC_CUTOFF"
  fi
fi

########################
### STEP 02 - Cicero ###
########################
echo "Step 02 - $(date +'%Y.%m.%d %H:%M:%S') - Cicero"
{
prepareCiceroInput.pl -o $CICERO_DATADIR/$SAMPLE -genome $GENOME -p $SAMPLE -l $LEN -s 250 -f $CICERO_DATADIR/$SAMPLE/${SAMPLE}.gene_info.txt $sc_cutoff_arg
for SCFILE in $CICERO_DATADIR/$SAMPLE/$SAMPLE.*.SC; do
    echo "Cicero.pl -i $BAMFILE -o $(echo $SCFILE | sed 's/\.SC//g') -l $LEN -genome $GENOME -f $SCFILE $cluster_arg"
done > cmds-02.sh
parallel --joblog 02_Cicero.log $PARALLEL_ARG < cmds-02.sh
} 1> 02_Cicero.out 2> 02_Cicero.err

#########################
### STEP 03 - Combine ###
#########################
echo "Step 03 - $(date +'%Y.%m.%d %H:%M:%S') - Combine"
{
cat $CICERO_DATADIR/$SAMPLE/*/unfiltered.fusion.txt | sort -V -k 9,9 -k 10,10n -k 11,11n > $CICERO_DATADIR/$SAMPLE/unfiltered.fusion.txt
cat $CICERO_DATADIR/$SAMPLE/*/unfiltered.internal.txt | sort -V -k 9,9 -k 10,10n -k 11,11n > $CICERO_DATADIR/$SAMPLE/unfiltered.internal.txt
} 1> 03_Combine.out 2> 03_Combine.err

## QC
if [ ! -s $CICERO_DATADIR/$SAMPLE/unfiltered.fusion.txt ]
then
  echo "No fusions in unfiltered data."
fi
if [ ! -s $CICERO_DATADIR/$SAMPLE/unfiltered.internal.txt ]
then
  echo "No internal events in unfiltered data."
fi 

##########################
### STEP 04 - Annotate ###
##########################
echo "Step 04 - $(date +'%Y.%m.%d %H:%M:%S') - Annotate"
{
echo "annotate.pl -i $BAMFILE -o $CICERO_DATADIR/$SAMPLE -genome $GENOME -l $LEN -s $SAMPLE -f $CICERO_DATADIR/$SAMPLE/${SAMPLE}.gene_info.txt -j $JUNCTIONS $cluster_arg" > cmds-04.sh
echo "annotate.pl -i $BAMFILE -o $CICERO_DATADIR/$SAMPLE -genome $GENOME -l $LEN -s $SAMPLE -f $CICERO_DATADIR/$SAMPLE/${SAMPLE}.gene_info.txt -internal $cluster_arg" >> cmds-04.sh
parallel --joblog 04_Annotate.log $PARALLEL_ARG < cmds-04.sh

mv $CICERO_DATADIR/$SAMPLE/blacklist.new.txt $CICERO_DATADIR/$SAMPLE/blacklist.new.fusions.txt
cat $CICERO_DATADIR/$SAMPLE/annotated.fusion.txt $CICERO_DATADIR/$SAMPLE/annotated.internal.txt > $CICERO_DATADIR/$SAMPLE/annotated.all.txt
} 1> 04_Annotate.out 2> 04_Annotate.err

## QC 
if [ $(wc -l $CICERO_DATADIR/$SAMPLE/annotated.all.txt | awk '{ print $1 }') -eq 1 ]
then
  echo "No annotated events found"
fi

########################
### STEP 05 - Filter ###
########################
echo "Step 05 - $(date +'%Y.%m.%d %H:%M:%S') - Filter"
{
cicero_filter.sh $CICERO_DATADIR $SAMPLE $GENOME
cp $CICERO_DATADIR/$SAMPLE/final_fusions.txt $CICERO_DATADIR/$SAMPLE/final_fusions.report.html
} 1> 05_Filter.out 2> 05_Filter.err

## QC
if  [ $(wc -l $CICERO_DATADIR/$SAMPLE/final_fusions.txt | awk '{ print $1 }') -eq 1 ]
then
  echo "No events found"
fi

############################
### Kill the blat server ###
############################
if [[ $BLAT_SERVER_PID ]]; then
    echo "Killing blat server with pid $BLAT_SERVER_PID" 2>&1
    kill $BLAT_SERVER_PID
fi

#############################
### Blow away tmp configs ###
#############################
#rm -rf $SJ_CONFIGS

