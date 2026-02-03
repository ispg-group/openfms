#!/bin/bash

# Must have MPI, TeraChem, and FMS90 in environment
# Must pass GPUS, TCEXE, and FMSEXE in on commandline as env variables
# Example: GPUS=0,1 TCEXE=terachem FMSEXE=FMS.tc.e ./run.sh
TCEXE=terachem
FMSEXE=openfms.tc
GPUS=0

MPI_ADDITIONAL_ARGUMENT="-nameserver localhost"

# Generate random port number to avoid conflicts
server=tcfms_port$(( ( RANDOM % 10000 ) + 1 ))
printf "&tc\nserver_name = '%s'\n/\n" "$server" > tc_input

if pgrep hydra_nameserve; then hydra_nameserver & sleep 1; fi

# For MPICH, OpenFabrics interface works in general
export MPIR_CVAR_CH4_NETMOD=ofi

# mpiexec.hydra seems to be the most effective way to start the two processes
mpiexec.hydra "$MPI_ADDITIONAL_ARGUMENT" -n 1 "$TCEXE" -g "$GPUS" -U2 --MPIPort="$server" &> tc.out &
PID_TC=$!
sleep 2 # grace time for terachem initialization (doesn't involve GPU initialization so should be fast)

mpiexec.hydra "$MPI_ADDITIONAL_ARGUMENT" -n 1 "$FMSEXE" &> fms.out &
PID_FMS=$!

echo "Both FMS(pid=$PID_FMS) and TC(pid=$PID_TC) started. Waiting for them to finish..."
echo "(Monitor tc.out and FMS.out for progress)"
# Should be replace with "wait -n" once we have bash 4.3
# Note about 'kill -0' https://unix.stackexchange.com/questions/169898/what-does-kill-0-do
while ( (kill -0 $PID1 >& /dev/null) && (kill -0 $PID2 >& /dev/null) ); do sleep 1; done
sleep 5 # grace time for program termination

# If one dies and the other doesn't, kill the other.
# This logic will be triggered if one dies before the other even starts,
# or if MPI doesn't close the other properly. Both cases are observed in practice.
# The following logic is not robust against OS pid reuse!
if ! (kill -0 $PID1 >& /dev/null); then
    wait $PID1; RETURN1=$?; echo "TC exited ($RETURN1)";
    if [ $RETURN1 -ne 0 ]; then
        kill $PID2; wait $PID2; RETURN2=$?; echo "FMS killed ($RETURN2)";
    fi
elif ! (kill -0 $PID2 >& /dev/null); then
    wait $PID2; RETURN2=$?; echo "FMS exited ($RETURN2)";
    if [ $RETURN2 -ne 0 ]; then
        kill $PID1; wait $PID1; RETURN1=$?; echo "TC killed ($RETURN1)";
    fi
fi

