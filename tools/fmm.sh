#!/usr/bin/env bash

FMM_MY_RANK=${PMIX_RANK}

# format the rank id
FMM_NCHARS_RANK=$(echo -n $((FMM_NTOTRANKS-1)) | wc -c)
FMM_RANK_ID_OUT=`printf "%0${FMM_NCHARS_RANK}d\n" "${FMM_MY_RANK}"`

# read the last argument
FMM_OUTPUT_FILE=${!#}

# read the optional argument #8
if [ $# -ge 9 ] ; then
  FMM_PRT_RANK=$8
else
  FMM_PRT_RANK=''
fi

# read the optional argument #9
if [ $# -ge 10 ] ; then
  FMM_SUBCOMSZ=$9
else
  FMM_SUBCOMSZ=''
fi

# read the optional argument #10
if [ $# -ge 11 ] ; then
  FMM_PLAN_RANKS=${10}
else
  FMM_PLAN_RANKS=''
fi

# read the optional argument #11
if [ $# -ge 12 ] ; then
  FMM_DD_LEVEL=${11}
else
  FMM_DD_LEVEL=''
fi

if [[ -z "${FMM_OPENMP_DEBUG}" ]]; then
  # this is a run without OpenMP debug
  unset LIBOMPTARGET_DEBUG

  if [ "${FMM_PRT_RANK}" -ge "0" ] ; then
    # FMM_PRT_RANK has a non-negative value

    # make sure that FMM_PRT_RANK doesn't exceed the numer of ranks
    if [ "${FMM_PRT_RANK}" -ge "${FMM_NTOTRANKS}" ] ; then
      FMM_PRT_RANK=$((FMM_NTOTRANKS-1))
    fi

    # only one MPI rank FMM_PRT_RANK writes to stdout
    ./fmm.x $1 $2 $3 $4 $5 $6 $7 ${FMM_PRT_RANK} ${FMM_SUBCOMSZ} ${FMM_PLAN_RANKS} ${FMM_DD_LEVEL} >> ${FMM_OUTPUT_FILE} 2>&1
  else
    # FMM_PRT_RANK has a negative value
    # this is an MPI debug mode
    # all MPI ranks write to stdout redirected into a separate output file
    ./fmm.x $1 $2 $3 $4 $5 $6 $7 ${FMM_PRT_RANK} ${FMM_SUBCOMSZ} ${FMM_PLAN_RANKS} ${FMM_DD_LEVEL} >> ${FMM_OUTPUT_FILE}.${FMM_RANK_ID_OUT} 2>&1
  fi
else
  # this is a run with OpenMP debug
  echo "FMM_OPENMP_DEBUG=${FMM_OPENMP_DEBUG}"

  # obtain backtraces from all MPI processes after 60 seconds of hanging
  export MPIR_CVAR_PROGRESS_TIMEOUT=60

  if [ "${FMM_PRT_RANK}" -ge "0" ] ; then
    # debug a single MPI rank when FMM_PRT_RANK has a non-negative value

    # make sure that FMM_PRT_RANK doesn't exceed the numer of ranks
    if [ "${FMM_PRT_RANK}" -ge "${FMM_NTOTRANKS}" ] ; then
      FMM_PRT_RANK=$((FMM_NTOTRANKS-1))
    fi

    if [ "${FMM_MY_RANK}" -eq "${FMM_PRT_RANK}" ] ; then
      # enable OpenMP debug for a single MPI rank FMM_PRT_RANK
      export LIBOMPTARGET_DEBUG=1
      ./fmm.x $1 $2 $3 $4 $5 $6 $7 ${FMM_PRT_RANK} ${FMM_SUBCOMSZ} ${FMM_PLAN_RANKS} ${FMM_DD_LEVEL} >> ${FMM_OUTPUT_FILE}.${FMM_RANK_ID_OUT} 2>&1
    else
      # do not save output from other MPI ranks
      ./fmm.x $1 $2 $3 $4 $5 $6 $7 ${FMM_PRT_RANK} ${FMM_SUBCOMSZ} ${FMM_PLAN_RANKS} ${FMM_DD_LEVEL} >> /dev/null 2>&1
    fi

  else
    # apply OpenMP debug to all MPI ranks
    export LIBOMPTARGET_DEBUG=1

    # write the debug output from each MPI rank into a separate file
    ./fmm.x $1 $2 $3 $4 $5 $6 $7 ${FMM_PRT_RANK} ${FMM_SUBCOMSZ} ${FMM_PLAN_RANKS} ${FMM_DD_LEVEL} >> ${FMM_OUTPUT_FILE}.${FMM_RANK_ID_OUT} 2>&1
  fi
fi
