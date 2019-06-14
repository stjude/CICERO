#!/usr/bin/env bash

CICERO_ROOT=`readlink -f $(dirname ${BASH_SOURCE[0]})/../../`
# TODO - need to totall revamp the config files
export SJ_CONFIGS=$CICERO_ROOT/configs.${HOSTNAME}.$$.tmp
cp -rf $CICERO_ROOT/configs $SJ_CONFIGS
echo "SJ_CONFIGS=$SJ_CONFIGS"

MAX_PORT_NUM=32000
MIN_PORT_NUM=2000

BAMFILE=""
JUNCTIONS=""
REFFA=""
GENOME=""
SAMPLE=""
TARGET="TRANSCRIPTOME"
NCORES=1
OUTDIR=$(pwd)

usage() {
    echo "Cicero [-h] [-n ncores] -b bamfile -g genome -r refdir [-j junctions] [-o outdir]"
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
    esac
    shift
done

#######################
### Validate inputs ###
#######################
if [[ ! -f $BAMFILE ]]; then
    >&2 echo "ERROR: Bamfile '$BAMFILE' DNE"
    usage
    exit 1
elif [[ ! -f $BAMFILE.bai ]]; then
    >&2 echo "ERROR: Bam index (.bai) '$BAMFILE.bai' DNE"
    usage
    exit 1
elif [[ ! -d $REFDIR ]]; then
    >&2 echo "ERROR: Reference directory '$REFFA' DNE"
    usage
    exit 1
elif [[ ! $GENOME ]]; then
    >&2 echo "ERROR: GENOME not defined"
    usage
    exit 1
elif [[ $JUNCTIONS && ! -f $JUNCTIONS ]]; then
    >&2 echo "ERROR: Junctions file '$JUNCTIONS' DNE"
    usage
    exit 1
fi

SAMPLE=$(basename $BAMFILE .bam)
if [[ $GENOME != "GRCh37-lite" ]]; then
    >&2 echo "ERROR: Only genome 'GRCh37-lite' is currently supported'"
    exit 1
fi

##############################
### Construct config files ###
##############################
for CONFIG_TYPE_DIR in $SJ_CONFIGS/*; do
    CONFIG_TYPE=$(basename $CONFIG_TYPE_DIR)
    for CONFIG_TEMPLATE in $CONFIG_TYPE_DIR/*.template.txt; do
        CONFIG_FILE=$CONFIG_TYPE_DIR/$(basename $CONFIG_TEMPLATE .template.txt).config.txt
        while read KEY VALUE; do
            if [[ $(echo $VALUE | grep '/') ]]; then
                echo -e "$KEY\t$REFDIR$VALUE"
            elif [[ $KEY == BLAT_PORT && $VALUE == 0 ]]; then
                echo -e "$KEY\t$$"
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
# Check if the port is open, if not, choose a random number until we find a free port
while [[ $(netstat -tulpn 2> /dev/null | grep LISTEN | grep ":$BLAT_PORT\b") ]]; do
    # TODO: Add warning for port change
    echo "WARN: Port $BLAT_PORT is currently in use. Trying new port..."
    # Due to the range of RANDOM, we need to narrow it down using some math
    RANGE=$(($MAX_PORT_NUM-$MIN_PORT_NUM+1))    # Find the total possible numbers we want
    RAND_PORT_NUM=$RANDOM                       # Assign out the $RANDOM var 
    let "RAND_PORT_NUM %= $RANGE"               # Mod the $RANDOM var using the range
    BLAT_PORT=$(($RAND_PORT_NUM+$MIN_PORT_NUM)) # Add the resulting number to the min to get it within the range
done

# Set inputs
mkdir -p $OUTDIR
cd $OUTDIR
CICERO_RUNDIR=CICERO_RUNDIR
CICERO_DATADIR=CICERO_DATADIR
CICERO_CONFIG=CICERO_CONFIG

##################
### Start BLAT ###
##################

# Start blat server
gfServer status $BLAT_HOST $BLAT_PORT 1> ${GFSERVER_LOG}.out 2> ${GFSERVER_LOG}.err
RETURN_CODE=$?
if [[ $RETURN_CODE != 0 ]]; then
    gfServer start $BLAT_HOST $BLAT_PORT "$TWOBIT" -log=${GFSERVER_LOG}.log 1>> ${GFSERVER_LOG}.out 2>> ${GFSERVER_LOG}.err  &
    BLAT_SERVER_PID=$!
    RETURN_CODE=1
    echo "Starting local blat server"
    NUM_CHECKS=0
    echo $BLAT_SERVER_PID
    while [[ $RETURN_CODE != 0 && $NUM_CHECKS -le 20 ]]; do
        echo "Blat server not yet running ..."
        if  ps -p $BLAT_SERVER_PID > /dev/null; then
            gfServer status $BLAT_HOST $BLAT_PORT 1>> ${GFSERVER_LOG}.out 2>> ${GFSERVER_LOG}.err
            RETURN_CODE=$?
            NUM_CHECKS=$(($NUM_CHECKS + 1))
            sleep 5
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
echo "Step 01 - $(date +'%Y.%m.%d %H:%M:%S') - ExtractSClips"
{
LEN=$(getReadLength.sh $BAMFILE)
get_sc_cmds.pl -i $BAMFILE -genome $GENOME -o $CICERO_DATADIR/$SAMPLE -l $LEN > cmds-01.sh
echo "get_geneInfo.pl -i $BAMFILE -o $CICERO_DATADIR/$SAMPLE -l $LEN -genome $GENOME -s $SAMPLE" >> cmds-01.sh
parallel -j $NCORES < cmds-01.sh
} 1> 01_ExtractSClips.out 2> 01_ExtractSClips.err

########################
### STEP 02 - Cicero ###
########################
echo "Step 02 - $(date +'%Y.%m.%d %H:%M:%S') - Cicero"
{
prepareCiceroInput.pl -o $CICERO_DATADIR/$SAMPLE -genome $GENOME -p $SAMPLE -l $LEN -s 250 -f $CICERO_DATADIR/$SAMPLE/${SAMPLE}.gene_info.txt
for SCFILE in $CICERO_DATADIR/$SAMPLE/$SAMPLE.*.SC; do
    echo "Cicero.pl -i $BAMFILE -o $(echo $SCFILE | sed 's/\.SC//g') -l $LEN -genome $GENOME -f $SCFILE"
done > cmds-02.sh
parallel -j $NCORES < cmds-02.sh
} 1> 02_Cicero.out 2> 02_Cicero.err

#########################
### STEP 03 - Combine ###
#########################
echo "Step 03 - $(date +'%Y.%m.%d %H:%M:%S') - Combine"
{
cat $CICERO_DATADIR/$SAMPLE/*/unfiltered.fusion.txt > $CICERO_DATADIR/$SAMPLE/unfiltered.fusion.txt
cat $CICERO_DATADIR/$SAMPLE/*/unfiltered.internal.txt > $CICERO_DATADIR/$SAMPLE/unfiltered.internal.txt
} 1> 03_Combine.out 2> 03_Combine.err

##########################
### STEP 04 - Annotate ###
##########################
echo "Step 04 - $(date +'%Y.%m.%d %H:%M:%S') - Annotate"
{
annotate.pl -i $BAMFILE -o $CICERO_DATADIR/$SAMPLE -genome $GENOME -l $LEN -s $SAMPLE -f $CICERO_DATADIR/$SAMPLE/${SAMPLE}.gene_info.txt -j $JUNCTIONS
annotate.pl -i $BAMFILE -o $CICERO_DATADIR/$SAMPLE -genome $GENOME -l $LEN -s $SAMPLE -f $CICERO_DATADIR/$SAMPLE/${SAMPLE}.gene_info.txt -internal # -j $JUNCTIONS
cat $CICERO_DATADIR/$SAMPLE/annotated.fusion.txt $CICERO_DATADIR/$SAMPLE/annotated.internal.txt > $CICERO_DATADIR/$SAMPLE/annotated.all.txt
} 1> 04_Annotate.out 2> 04_Annotate.err

########################
### STEP 05 - Filter ###
########################
echo "Step 05 - $(date +'%Y.%m.%d %H:%M:%S') - Filter"
{
cicero_filter.sh $CICERO_DATADIR $SAMPLE $GENOME
cp $CICERO_DATADIR/$SAMPLE/final_fusions.txt $CICERO_DATADIR/$SAMPLE/final_fusions.report.html
} 1> 05_Filter.out 2> 05_Filter.err

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
rm -rf $SJ_CONFIGS

