#!/bin/bash

# Generate cubic i,j,k indices based on provided number of replication layers

if [ $# != 6 ] ; then
  echo ""
  echo "Usage: $0 nRepLayers datasetID selectionID xStart yStart zStart"
  exit
fi

# Number of replication layers
NREPLAYERS=$1
NREPLAYERS1=$((NREPLAYERS + 1))

# Number of replicas
NREPLICAS=$((NREPLAYERS1 * NREPLAYERS1 * NREPLAYERS1))

# dataset id
DATASETID=$2

# selection id
SELECTIONID=$3

# starting values of the replication index
XSTART=$4
YSTART=$5
ZSTART=$6

#-------------------------------------------------
echo "Number of replication layers = ${NREPLAYERS}"
echo "Number of replicas           = ${NREPLICAS}"
echo "Number of selected groups    = ${SELECTIONID}"
echo ""

# nDataSets, nReplicas, nGroups
echo "1 ${NREPLICAS} ${SELECTIONID}"
echo "datasetFilename"

Z=${ZSTART}
for(( k=0; k<=NREPLAYERS; k++)); do
  Y=${YSTART}
  for(( j=0; j<=NREPLAYERS; j++)); do
    X=${XSTART}
    for(( i=0; i<=NREPLAYERS; i++)); do
       echo "$DATASETID $SELECTIONID   $X $Y $Z"    # dataset_id, selection_group,  i, j, j
       ((X++))
    done
    ((Y++))
  done
  ((Z++))
done
