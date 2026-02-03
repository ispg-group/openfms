#!/bin/bash

# This is a script for launching an AIMS simulations with TeraChem MPI interface.
#
# You must ensure that the necessary TeraChem environment is setup here,
# a path to OpenFMS binary, and importantly MPI environment that was used to compile both of them.

TCEXE=terachem
FMSEXE=openfms.tc
MPIRUN=mpiexec.hydra
GPUS=0

MPI_ADDITIONAL_ARGUMENTS='-nameserver localhost'

### Let's check all the binaries are available first
if ! which "$MPIRUN" > /dev/null; then
  echo "ERROR: Executable $MPIRUN not found"
  exit 1
fi
if ! which "$TCEXE" > /dev/null; then
  echo "ERROR: Executable $TCEXE not found"
  exit 1
fi
if ! which "$FMSEXE" > /dev/null; then
  echo "ERROR: Executable $FMSEXE not found"
  exit 1
fi

# Generate random port number to avoid conflicts
server=tcfms_port$(( ( RANDOM % 10000 ) + 1 ))
printf "&tc\nserver_name = '%s'\n/\n" "$server" > tc_input

if ! pgrep -f hydra_nameserver; then hydra_nameserver & sleep 1; fi

# For MPICH, OpenFabrics interface works in general
export MPIR_CVAR_CH4_NETMOD=ofi

# shellcheck disable=SC2086
$MPIRUN $MPI_ADDITIONAL_ARGUMENTS -n 1 "$TCEXE" -g "$GPUS" -U2 --MPIPort="$server" &> tc.out &
PID_TC=$!
sleep 2 # grace time for terachem initialization (doesn't involve GPU initialization so should be fast)

# shellcheck disable=SC2086
$MPIRUN $MPI_ADDITIONAL_ARGUMENTS -n 1 "$FMSEXE" &> fms.out &
PID_FMS=$!

echo "Both OpenFMS(pid=$PID_FMS) and TeraChem(pid=$PID_TC) have started, waiting for them to finish..."
echo "(Monitor tc.out, fms.out and FMS.out for progress)"
# Should be replace with "wait -n" once we have bash 4.3
# Note about 'kill -0' https://unix.stackexchange.com/questions/169898/what-does-kill-0-do
while ( (kill -0 $PID_TC >& /dev/null) && (kill -0 $PID_FMS >& /dev/null) ); do sleep 1; done
sleep 5 # grace time for program termination

# If one dies and the other doesn't, kill the other.
# This logic will be triggered if one dies before the other even starts,
# or if MPI doesn't close the other properly. Both cases are observed in practice.
# The following logic is not robust against OS pid reuse!
if ! (kill -0 $PID_TC >& /dev/null); then
    wait $PID_TC; RETURN1=$?; echo "TC exited ($RETURN1)";
    if [ $RETURN1 -ne 0 ]; then
        kill $PID_FMS; wait $PID_FMS; RETURN2=$?; echo "FMS killed ($RETURN2)";
    fi
elif ! (kill -0 $PID_FMS >& /dev/null); then
    wait $PID_FMS; RETURN2=$?; echo "FMS exited ($RETURN2)";
    if [ $RETURN2 -ne 0 ]; then
        kill $PID_TC; wait $PID_TC; RETURN1=$?; echo "TC killed ($RETURN1)";
    fi
fi
