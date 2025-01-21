#!/bin/bash

#source setup_env.sh
source /gapps/root/Linux_RHEL7-x86_64-gcc4.8.2/root-6.18.00/bin/thisroot.csh
# only root needed

RUNNUM=${1-none}
#MAXEVT=${2-0}

if [[ ${RUNNUM} == "none" ]] ; then
    echo "================================="
    echo " Usage: ./$0 <RunNum> [Max_Events] "
    echo "================================="
    exit 0;
fi

echo "====>  Process RUN=$RUNNUM <=========="

root --web=off -l <<EOC
.x trd_mlp_halld24.C++(${RUNNUM})
EOC
