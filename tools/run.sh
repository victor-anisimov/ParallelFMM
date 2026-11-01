#!/bin/bash

NARGS_REQUIRED=17

if [ $# -ne ${NARGS_REQUIRED} ] ; then
  echo "Error: Wrong number of arguments $# supplied to ${0};   Required ${NARGS_REQUIRED} arguments"
  exit 1
else

  # job parameters
  export FMM_NUMNODES=${1}
  export FMM_WALLTIME=${2}
  export FMM_RUNQUEUE=${3}
  export FMM_JOB_NAME=${4}
  export FMM_NODE_RANKS=${5}

  # program arguments
  export FMM_INPUT=${6}
  export FMM_PRT_RANK=${7}
  export FMM_SUBCOMSZ=${8}
  export FMM_PROBSIZE=${9}
  export FMM_SIDE_LEN=${10}
  export FMM_EXPAN_SZ=${11}
  export FMM_WS=${12}

  # batch arguments
  export FMM_RUN_START=${13}
  export FMM_RUN_END=${14}
  export FMM_RUN_INCREMENT=${15}

  # FMM plan queue
  export FMM_PLAN_RANKS=${16}

  # FMM ddLevel
  export FMM_DD_LEVEL=${17}

  if [ ${FMM_EXPAN_SZ} == "scan" ] ; then
    FMM_EXPAN_SZ="scan from ${FMM_RUN_START} to ${FMM_RUN_END} with increment ${FMM_RUN_INCREMENT}"
    FMM_SCAN_EXPAN_SZ=1
  else
    FMM_SCAN_EXPAN_SZ=0
  fi

  echo "Number of nodes          ${FMM_NUMNODES}"
  echo "Walltime                 ${FMM_WALLTIME}"
  echo "Queue                    ${FMM_RUNQUEUE}"
if [ -n "$FMM_ALLOCATION" ] ; then
  echo "Allocation               $FMM_ALLOCATION"
  export FMM_ALLOCATION="#PBS -A ${FMM_ALLOCATION}"
else
  echo "Allocation               undefined"
  exit
fi
if [ -n "$FMM_FILESYSTEM" ] ; then
  echo "Filesystems              $FMM_FILESYSTEM"
  export FMM_FILESYSTEM="#PBS -l filesystems=${FMM_FILESYSTEM}"
fi
  echo "Output job prefix        ${FMM_JOB_NAME}"
  echo "Number of ranks per node ${FMM_NODE_RANKS}"
  echo "FMM Plan nranks per node ${FMM_PLAN_RANKS}"
  echo "FMM ddLevel input value  ${FMM_DD_LEVEL}"
  echo "Input file name          ${FMM_INPUT}"
  echo "MPI print rank           ${FMM_PRT_RANK}"
  echo "Sub-communicator size    ${FMM_SUBCOMSZ}"
  echo "Problem size             ${FMM_PROBSIZE}"
  echo "Cube side-length         ${FMM_SIDE_LEN}"
  echo "Expansion size           ${FMM_EXPAN_SZ}"
  echo "Well-separation          ${FMM_WS}"
  echo "Run id start number      ${FMM_RUN_START}"
  echo "Run id end number        ${FMM_RUN_END}"
  echo "Run id increment         ${FMM_RUN_INCREMENT}"

  # determine the locality of node placement
  if [ ${FMM_NUMNODES} -le 8 ] ; then
    # place nodes on a chassis (8 nodes)
    FMM_NODE_PLACEMENT="-l place=group=tier1"  # one can also use "switch" instead of "tier1"
  else
    # place nodes on a rack (64 nodes) if the job fits in a single rack
    FMM_NODE_PLACEMENT="-l place=group=tier0"  # will be ignored if the number of nodes > 64
  fi

  for (( i=${FMM_RUN_START}; i<=${FMM_RUN_END}; i+=${FMM_RUN_INCREMENT} )) ; do

    # scan the expansion size
    if [ ${FMM_SCAN_EXPAN_SZ} -eq 1 ] ; then
      export FMM_EXPAN_SZ=${i}
    fi

    # formatted parameters
    export FMM_NODE_RANKS_OUT=`printf '%02d\n' "${FMM_NODE_RANKS}"`
    export FMM_NUMNODES_OUT=`printf '%04d\n' "${FMM_NUMNODES}"`
    export FMM_PROBSIZE_OUT=`printf '%03d\n' "${FMM_PROBSIZE}"`
    export FMM_SUBCOMSZ_OUT=`printf '%03d\n' "${FMM_SUBCOMSZ}"`
    export FMM_CU_OUT=`printf '%04.1f\n' "${FMM_SIDE_LEN}"`
    export FMM_NL_OUT=`printf '%02d\n' "${FMM_EXPAN_SZ}"`

    # formatted job id
    export FMM_JOB_ID=`printf '%02d\n' "$i"`
    # PBS job script name
    FMM_RUN_PBS="run.pbs.${FMM_JOB_ID}"
    # generate PBS job script
    printf '%s\n\n%s\n%s%s\n' "#!/bin/bash" "#PBS -l select=${FMM_NUMNODES} ${FMM_NODE_PLACEMENT}" \
           "$(cat tools/template.pbs)" ".${FMM_JOB_ID}" > ${FMM_RUN_PBS}
    # substitute parameters with the actual values
    envsubst '${FMM_PRT_RANK} ${FMM_SUBCOMSZ}     ${FMM_SUBCOMSZ_OUT} ${FMM_NODE_RANKS} \
              ${FMM_PROBSIZE} ${FMM_PROBSIZE_OUT} ${FMM_INPUT}        ${FMM_NODE_RANKS_OUT} \
              ${FMM_SIDE_LEN} ${FMM_EXPAN_SZ}     ${FMM_WS}           ${FMM_NL_OUT} \
              ${FMM_WALLTIME} ${FMM_RUNQUEUE}     ${FMM_JOB_NAME}     ${FMM_CU_OUT} \
              ${FMM_JOB_ID}   ${FMM_NUMNODES_OUT} ${FMM_OPENMP_DEBUG} ${FMM_PLAN_RANKS} \
	      ${FMM_FILESYSTEM} ${FMM_ALLOCATION} ${FMM_DD_LEVEL}' < ${FMM_RUN_PBS} > .tmp
    # rename the job script
    mv .tmp ${FMM_RUN_PBS}
    # print the job id to stdout
    echo ${FMM_JOB_ID}
    # submit the PBS job
    qsub ${FMM_RUN_PBS}
    # increment the value of expansion size
    #export FMM_EXPAN_SZ=$((FMM_EXPAN_SZ+1))
   done

fi
