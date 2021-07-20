#!/usr/bin/env bash

set -ex

VERSION=$1
if [ -z "$VERSION" ]
then
  VERSION="latest"
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "Script dir: $SCRIPT_DIR"

echo "Using CICERO version: $VERSION"

GENOME=GRCh38_no_alt
BAM=data/input/test.bam
REFDIR=reference
JUNCTIONS=data/input/test.bam.junctions.tab.shifted.tab

echo "Genome: $GENOME"
echo "BAM: $BAM"
echo "Reference directory: $REFDIR"
echo "Junctions file: $JUNCTIONS"

if [ ! -d $SCRIPT_DIR/$REFDIR ]
then
  echo "Fetching reference files"
  curl -O https://zenodo.org/record/5088371/files/reference.tar.gz
  echo "c2a35c360539c18c0f42cb26c686cb0e  reference.tar.gz" > reference.md5
  md5sum -c reference.md5
  tar -zxf reference.tar.gz
  find reference
else
  echo "Reference exists"
fi

echo "docker run -v ${SCRIPT_DIR}/${REFDIR}:/reference -v ${SCRIPT_DIR}:/data ghcr.io/stjude/cicero:${VERSION} Cicero.sh -b /data/${BAM} -g ${GENOME} -r /reference -o /data/output -j /data/${JUNCTIONS}"
docker run --memory="2.5g" --memory-swap="10g" -v ${SCRIPT_DIR}/${REFDIR}:/reference -v ${SCRIPT_DIR}:/data ghcr.io/stjude/cicero:${VERSION} Cicero.sh -b /data/${BAM} -g ${GENOME} -r /reference -o /data/output -j /data/${JUNCTIONS}
#docker run -v ${SCRIPT_DIR}/${REFDIR}:/reference -v ${SCRIPT_DIR}:/data ghcr.io/stjude/cicero:${VERSION} Cicero.sh -b /data/${BAM} -g ${GENOME} -r /reference -o /data/output -j /data/${JUNCTIONS}

