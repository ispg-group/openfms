#!/bin/bash

set -u

rm -f errors.dat err.out FMS.out fms.out
if [[ $1 = "clean" ]];then
  exit 0
fi

FMSEXE=$1

# Test that ABIN fails with invalid cmdline argument
$FMSEXE -invalid > fms.out 2>> err.out

# Test that ABIN prints help, without an error
$FMSEXE -h >> fms.out 2>> err.out || echo "ERROR when printing help" >> error.out
$FMSEXE --help >> fms.out 2>> err.out || echo "ERROR when printing help" >> error.out

# Test that ABIN prints version, without an error
$FMSEXE -v >> fms.out 2>> err.out || echo "ERROR when printing version" >> error.out
$FMSEXE --version >> fms.out 2>> err.out || echo "ERROR when printing version" >> error.out

# Filter out "STOP X" or "X" from stderr, where X is the return code
# This hack is tested for gfortran and intel compilers
grep -E -v '^(STOP )?[01]$' err.out > errors.dat
