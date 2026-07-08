#!/bin/bash

# Stop on any error
set -euo pipefail

# Make sure that '[A-Z]' regex works as expected,
# i.e. only matches capital letters, regardless of global collate settings.
# This workaround is needed for MacOS 15,
# but is probably a good idea in general to make this script robust.
export LC_ALL=C

# Parameters are passed from Makefile.
NUMPARAM=3
if [[ $# -ne $NUMPARAM ]]; then
  echo "ERROR: Incorrect number of parameters passed to $0"
  echo "Invoked as:"
  echo "$0 $@"
  exit 1
fi

if [[ "$1" = /* ]]; then
  FMSEXE="$1"
else
  FMSEXE="$PWD/$1"
fi
FMSOUT="err.out"
TESTS="$2"
ACTION="$3"

if [[ ! -f $FMSEXE ]]; then
    echo "ERROR: File $FMSEXE does not exist!"
    exit 2
fi

if [[ $ACTION != "makeref" && $ACTION != "clean" && $ACTION != 'test' ]]; then
    echo "ERROR: Parameter '$ACTION' must be one of 'test', 'makeref' or 'clean'"
    exit 2
fi

if [[ $ACTION = "makeref" && $TESTS = "all" ]];then
  echo "ERROR: You should not call makeref on all tests at once."
  echo "Specify a concrete test which you want to modify, e.g."
  echo "make makeref TEST=CMD"
  exit 1
fi

cd $(dirname $0)
TESTDIR=$PWD

function diff_files {
  local return_status=0
  local error_code
  local ref_file
  local test_file
  # Compare test results with all existing reference files
  local reference_files=$(ls *.ref)
  if [[ -z $reference_files ]];then
    echo "ERROR: No reference files were found"
    return 1
  fi

  for ref_file in $reference_files
  do
    test_file=$(basename $ref_file .ref)
    if [[ ! -f $test_file ]];then
      # The output file does not exist.
      # Something went seriously wrong, we probably crashed prematurely.
      # No need for further checks, exit NOW.
      echo "ERROR: Could not find output file \"$test_file\""
      return 1
    fi

    error_code=0
    diff -q $test_file $ref_file > /dev/null || error_code=$?
    if [[ $error_code -ne 0 ]];then

       if [[ $test_file = "errors.dat" ]]; then
          echo "Files $ref_file and $test_file differ!"
          diff --color=always $ref_file $test_file | tee $test_file.diff
          return_status=$error_code
          continue
       fi

       # The reference file is different, but maybe it's just numerical noise?
       error_code=0
       diff -y -W 500  $test_file $ref_file | grep -e '|' -e '<' -e '>' > $test_file.diff

       ../numdiff.py $test_file.diff || error_code=$?

       if [[ $error_code -ne 0 ]];then
          # The changes were bigger that the thresholds specified in numdiff.py
          return_status=$error_code
       fi
    fi
  done
  return $return_status
}

# Yes, I know this is silly, let me have my OCD okay?
function num_tests {
  local count=$1
  if [[ count -eq 1 ]]; then
    echo "1 test"
  else
    echo "$count tests"
  fi
}

# Update already existing reference files.
# Called by `make makeref TEST=TEST_FOLDER`
# If you're creating a completely new test,
# you need to create the reference files manually.
function makeref {
  local ref_file
  local test_file
  local reference_files=$(ls *.ref)
  echo "Making new reference files."
  if [[ -z $reference_files ]];then
    echo "ERROR: No reference files were found."
    exit 1
  fi
  for ref_file in $reference_files
  do
    test_file=$(basename $ref_file .ref)
    if [[ ! -f $test_file ]];then
      echo "ERROR: Could not find output file \"$test_file\""
      exit 1
    fi
    echo "mv $test_file $ref_file"
    mv $test_file $ref_file
  done
}

function clean {
  if [[ -f "test.sh" ]];then
    ./test.sh clean
  fi
  rm -rf $*
  rm -f *.diff
}

# List of all possible FMS output files.
# Used by `make testclean` to cleanup test directories.
# TODO!
output_files=(?.dat CFxn.dat SDot.dat fort.0 Out.xyz forces.*.xyz positions.*.xyz positions.*.dcd *.[1-9] *.out *.log Checkpoint.txt Last_Bundle.txt)

# Run all tests
if [[ $TESTS = "all" ]];then

   folders=$(ls -d [A-Z]*)

else

   # Only one test selected, e.g. by running
   # make test TEST=CMD
   folders=${TESTS}

fi

echo "Running tests in directories:"
echo ${folders}

errors=0
passed=0

for dir in ${folders[@]}
do
   if [[ ! -d $dir ]];then
      echo "Directory $dir not found. Exiting prematurely."
      exit 1
   fi
   echo -en "Running $dir\t"
   cd $dir

   # Always clean the test directory before runnning the test.
   clean ${output_files[@]}

   # If we just want to clean the directories,
   # we skip the the actual test here
   if [[ $ACTION = "clean" ]];then
      echo "Cleaning files in directory $dir"
      cd $TESTDIR
      continue
   fi

   # For tests that require restart
   if [[ -f Last_Bundle.txt.initial ]]; then
      cp Last_Bundle.txt.initial Last_Bundle.txt
   fi

   # Don't exit on errors in the rest of the script, since fms invocation can fail
   # and this makes the code simpler.
   set +e
   if [[ -f test.sh ]]; then
     # Run custom testing script
     ./test.sh $FMSEXE
   else
     $FMSEXE &> $FMSOUT
   fi

   if [[ $ACTION = "makeref" ]];then

      makeref

   else

      if diff_files; then
        echo -e "\033[0;32mPASSED\033[0m"
        let passed++
      else
        let errors++
        # Uncomment these for debugging on GitHub
        # echo "=== FMS STDOUT & STDERR ==="
        # cat $FMSOUT
        echo -e "$dir \033[0;31mFAILED\033[0m"
      fi
   fi

   echo "---------------------------------------"

   cd $TESTDIR
done

echo " "

echo -e "\033[0;32m$(num_tests $passed) PASSED.\033[0m"
if [[ ${errors} -ne 0 ]];then
  echo -e "\033[0;31m$(num_tests $errors) FAILED\033[0m."
  exit 1
fi
