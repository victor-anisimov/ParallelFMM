#!/bin/bash

# job parameters
export FMM_NUMNODES=4
export FMM_WALLTIME=00:05:00
export FMM_RUNQUEUE=[job queue]
export FMM_ALLOCATION=[allocation project]
export FMM_FILESYSTEM=home:flare
export FMM_JOB_NAME=job
export FMM_NODE_RANKS=12

# FMM Plan queue
export FMM_PLAN_RANKS=${FMM_NODE_RANKS}
#export FMM_PLAN_RANKS=6

# program arguments
export FMM_INPUT=data/aon320

# setting FMM_PRT_RANK to a negative value enables stdout from each MPI rank into a separate file
export FMM_PRT_RANK=0
export FMM_SUBCOMSZ=1
export FMM_PROBSIZE=3
export FMM_SIDE_LEN=20
#export FMM_EXPAN_SZ=scan
export FMM_EXPAN_SZ=30
export FMM_WS=3
export FMM_DD_LEVEL=0

# batch arguments
export FMM_RUN_START=0
export FMM_RUN_END=0
export FMM_RUN_INCREMENT=1

# disable OpenMP debug, default
unset FMM_OPENMP_DEBUG

# enable OpenMP debug
#export FMM_OPENMP_DEBUG=1

./tools/run.sh         \
  ${FMM_NUMNODES}      \
  ${FMM_WALLTIME}      \
  ${FMM_RUNQUEUE}      \
  ${FMM_JOB_NAME}      \
  ${FMM_NODE_RANKS}    \
  ${FMM_INPUT}         \
  ${FMM_PRT_RANK}      \
  ${FMM_SUBCOMSZ}      \
  ${FMM_PROBSIZE}      \
  ${FMM_SIDE_LEN}      \
  ${FMM_EXPAN_SZ}      \
  ${FMM_WS}            \
  ${FMM_RUN_START}     \
  ${FMM_RUN_END}       \
  ${FMM_RUN_INCREMENT} \
  ${FMM_PLAN_RANKS}    \
  ${FMM_DD_LEVEL}
