#!/usr/bin/env bash

CICERO_ROOT=`readlink -f $(dirname ${BASH_SOURCE[0]})/../../`
# TODO - need to totall revamp the config files
export SJ_CONFIGS=$CICERO_ROOT/configs

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
for CONFIG_TYPE_DIR in $CICERO_ROOT/configs/*; do
    CONFIG_TYPE=$(basename $CONFIG_TYPE_DIR)
    for CONFIG_TEMPLATE in $CONFIG_TYPE_DIR/*.template.txt; do
        CONFIG_FILE=$CONFIG_TYPE_DIR/$(basename $CONFIG_TEMPLATE .template.txt).config.txt
        while read KEY VALUE; do
            if [[ $(echo $VALUE | grep '/') ]]; then
                echo -e "$KEY\t$REFDIR$VALUE"
            else
                echo -e "$KEY\t$VALUE"
            fi
        done < $CONFIG_TEMPLATE > $CONFIG_FILE
    done
done

FASTA=$(awk -F$'\t' '$1=="FASTA"{print $2}' $CICERO_ROOT/configs/genome/$GENOME.config.txt)
BLAT_HOST=$(awk -F$'\t' '$1=="BLAT_HOST"{print $2}' $CICERO_ROOT/configs/genome/$GENOME.config.txt)
BLAT_PORT=$(awk -F$'\t' '$1=="BLAT_PORT"{print $2}' $CICERO_ROOT/configs/genome/$GENOME.config.txt)
TWOBIT=$(awk -F$'\t' '$1=="TWOBIT"{print $2}' $CICERO_ROOT/configs/genome/$GENOME.config.txt)

##################
### Start BLAT ###
##################

# Start blat server
gfServer status $BLAT_HOST $BLAT_PORT > /dev/null 2>&1
RETURN_CODE=$?
if [[ $RETURN_CODE != 0 ]]; then
    gfServer start $BLAT_HOST $BLAT_PORT "$TWOBIT" -log=gfServer.log > /dev/null &
    BLAT_SERVER_PID=$!
    RETURN_CODE=1
    echo "Starting local blat server"
    while [[ $RETURN_CODE != 0 ]]; do
        echo "Blat server not yet running ..."
        gfServer status $BLAT_HOST $BLAT_PORT > /dev/null 2>&1
        RETURN_CODE=$?
        sleep 5
    done
    echo "Blat server is up!"
else
    echo "Blat server is already up!"
fi

# Set inputs
mkdir -p $OUTDIR
cd $OUTDIR
CICERO_RUNDIR=CICERO_RUNDIR
CICERO_DATADIR=CICERO_DATADIR
CICERO_CONFIG=CICERO_CONFIG

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
prepareCiceroInput.pl -o $CICERO_DATADIR/$SAMPLE -genome $GENOME -p $SAMPLE -l $LEN -s 250 -f $CICERO_DATADIR/$SAMPLE/$SAMPLE.gene_info.txt
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
annotate.pl -i $BAMFILE -o $CICERO_DATADIR/$SAMPLE -genome $GENOME -l $LEN -s $SAMPLE -f $CICERO_DATADIR/$SAMPLE/$SAMPLE.gene_info.txt -j $JUNCTIONS
annotate.pl -i $BAMFILE -o $CICERO_DATADIR/$SAMPLE -genome $GENOME -l $LEN -s $SAMPLE -f $CICERO_DATADIR/$SAMPLE/$SAMPLE.gene_info.txt -internal # -j $JUNCTIONS
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
    kill $BLAT_SERVER_PID
fi

