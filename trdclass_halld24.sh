#!/bin/bash


#source setup_env.sh
# source /gapps/root/Linux_RHEL7-x86_64-gcc4.8.2/root-6.18.00/bin/thisroot.sh
# only root needed

RUNNUM=${1-none}
MAXEVT=${2-0}
FRSTEVT=${3-0}
FILENUM=${4--1}

if [[ ${RUNNUM} == "none" ]] ; then
    echo "================================="
    echo " Usage: ./$0 <RunNum> [Max_Events] [First_Event] [File_Number]"
    echo "================================="
    exit 0;
fi

echo "=========>  Process RUN=$RUNNUM FILE=$FILENUM <=========="

root -l <<EOC
.L trdclass_halld24.C+g
trdclass_halld24 t(${RUNNUM},${MAXEVT},${FRSTEVT},${FILENUM})
t.Loop()
EOC
