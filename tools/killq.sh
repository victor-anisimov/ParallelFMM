#!/bin/bash

# $1 is PBS Jobname
# if $1 == NULL then kill all jobs
# if $1 != NULL then kill jobs that have Jobname == $1

for i in $(qstat -u ${USER} -f | awk -v val=$1 -v username=${USER} '{\
	if(val == NULL) \
	  {if($2 == username && $10 == "Q") print $0} \
	else \
	  {if($2 == username && $10 == "Q" && $4 == val) print $0} \
	  }' | awk -F. '{print $1}') ;
do
  qdel $i
  echo "$i killed"
done
