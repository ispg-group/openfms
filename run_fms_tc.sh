#!/bin/bash

# Must have intel, intel MPI, TeraChem, and FMS90 in environment
# Must pass GPUS, TCEXE, and FMSEXE in on commandline as env variables
# Example: GPUS=0,1 TCEXE=terachem FMSEXE=FMS.tc.e ./run.sh
echo "GPUS: $GPUS"
echo "TCEXE: $TCEXE"
echo "FMSEXE: $FMSEXE"
echo "MPI_ADDITIONAL_ARGUMENT: $MPI_ADDITIONAL_ARGUMENT (must be -nameserver localhost for intel MPI)"

# Generate random port number to avoid conflicts
server=tcfms_port$[ ( $RANDOM % 10000 ) + 1]
printf "&tc\nserver_name = '$server'\n" > tc_input

if [ -z $(ps aux | grep hydra_nameserver | grep -v "grep" | awk '{print $2}') ]; then hydra_nameserver & sleep 1; fi
# For MPICH, OpenFabrics interface works in general
export MPIR_CVAR_CH4_NETMOD=ofi
# mpiexec.hydra seems to be the most effective way to start the two processes
mpiexec.hydra $MPI_ADDITIONAL_ARGUMENT -n 1 $TCEXE -g $GPUS -U2 --MPIPort=$server > tc.out &
PID1=$!
sleep 5 # grace time for terachem initialization (FMS mode initialization doesn't include GPU startup, so it should be fast)
mpiexec.hydra $MPI_ADDITIONAL_ARGUMENT -n 1 $FMSEXE > fms.out &
PID2=$!

# Should be replace with "wait -n" once we have bash 4.3
while ( (kill -0 $PID1 >& /dev/null) && (kill -0 $PID2 >& /dev/null) ); do sleep 1; echo "Waiting..."; done
sleep 5 # grace time for program termination
# If one dies and the other doesn't, kill the other. This logic will be triggered if one dies before the other even starts, or if MPI doesn't close the other properly. Both cases are observed in practice.
# The following logic is not robust against OS pid resue!
if ! (kill -0 $PID1 >& /dev/null); then
    wait $PID1; RETURN1=$?; echo "1 exit, return = $RETURN1";
    if [ $RETURN1 -ne 0 ]; then
        kill $PID2; wait $PID2; RETURN2=$?; echo "2 killed, return = $RETURN2";
    fi
elif ! (kill -0 $PID2 >& /dev/null); then
    wait $PID2; RETURN2=$?; echo "2 exit, return = $RETURN2";
    if [ $RETURN2 -ne 0 ]; then
        kill $PID1; wait $PID1; RETURN1=$?; echo "1 killed, return = $RETURN1";
    fi
fi

