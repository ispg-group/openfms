#!/bin/bash

set -u

if [[ "$1" = "clean" ]]; then
  rm -f FMS.out err.out
  exit 0
fi

FMSEXE=$1
if [[ ! -f $FMSEXE ]]; then
  echo "ERROR when running test in $PWD"
  echo "FMS executable '$FMSEXE' does not exist!"
  exit 1
fi

# We only care about the error message here.
# Filtering out the other output helps to ensure
# consistency between compilers.
$FMSEXE 2>&1 | grep ERROR > err.out
