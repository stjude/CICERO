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

# Set inputs
mkdir -p $OUTDIR
cd $OUTDIR

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

