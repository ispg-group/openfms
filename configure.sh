#!/bin/bash
# Embarassigly barebones configure.sh, which generates the CONFIGFMS file used by make
# Takes in command-line options that define flags used during compilation and linking
# The first parameter to this script determines the interface that is compiled.
# Currently, only the MPI-based interface to TeraChem is available.
# By default, we only compile with the model potentials ('zero').


function show_help() {
echo "
NAME
   configure.sh  
DESCRIPTION
   Generates local compilation config file for OpenFMS and saves it to CONFIGFMS

USAGE
   configure [--interface] [--stdlib-lib path/to/stdlib-lib] [--stdlib-module path/to/stdlib-module]
OPTIONS
   --interface      Define electronic-structure interface with which OpenFMS is compiled.
                    Valid options:
                      'zero' for no interface, only model systems will be available. (default)
                      'tc' for TeraChem
   --stdlib-lib     location of STDLIB libraries (https://stdlib.fortran-lang.org/)
   --stdlib-module  location of STDLIB modules (https://stdlib.fortran-lang.org/)
   -h|--help        print this help and exit
"
}

function check_path() {
    if [[ ! -d $1 ]]; then
        echo "'$1' is not a directory!"
        exit 1
    fi
}


# Output file name
VARS=CONFIGFMS

ESP=zero

# Parse command-line options
PARAMS="$@"
while [[ $# -gt 0 ]]; do
    case "$1" in
      --interface)
        ESP=${2}
        if [[ $ESP != "tc" && $ESP != "zero" && $ESP != "quantics" ]]; then
            echo "ERROR: Unsupported interface target '$ESP'"
            echo "Supported targets: 'zero', 'tc'"
            exit 1
        fi
        shift 2
        ;;
      --stdlib-lib)
        STDLIB_LIBDIR=$2
        check_path "$STDLIB_LIBDIR"
        STDLIB=TRUE
        shift 2
        ;;
      --stdlib-module)
        STDLIB_MODULEDIR=$2
        if [[ $STDLIB_LIBDIR != "none" && $STDLIB_MODULEDIR = "none" ]]; then
            echo "ERROR: Please provide location of the STDLIB module directory"
            exit 1
        fi
        check_path "$STDLIB_MODULEDIR"
        if [[ $STDLIB_LIBDIR != "none" && $STDLIB_MODULEDIR != "none" ]]; then
            STDLIB=TRUE
        fi
        STDLIB=TRUE
        shift 2
        ;;
      --help|-h)
        show_help
        exit 0
        ;;
      *)
        echo "ERROR: Invalid option '$1'"
        show_help
        exit 1
        ;;
    esac
done

# Sanity check for command-line set variables
if [[ -z ${STDLIB_LIBDIR-} && ${STDLIB-} = "TRUE" ]]; then
    echo "ERROR: Please provide location of the STDLIB library directory"
    exit 1
fi
if [[ -z ${STDLIB_MODULEDIR-} && ${STDLIB-} = "TRUE" ]]; then
    echo "ERROR: Please provide location of the STDLIB module directory"
    exit 1
fi


if [[ -z $FC ]]; then
    if [[ $ESP = 'tc' ]];then
	# We need MPI compiler for TC interface
	# Look for intelmpi first, then gnu
	if [[ -n $(which mpiifort 2>/dev/null) ]]; then
	    FC=mpiifort
	elif [[ -n $(which mpifort 2>/dev/null) ]]; then
	    FC=mpifort
	else
	    echo "ERROR: Could not find MPI compiler for TC interface."
	    echo "(tried 'mpifort' and 'mpiifort')"
	    exit 1
	fi
    else
	if [[ -n $(which gfortran 2>/dev/null) ]]; then
	    FC=gfortran
	else
	    echo "ERROR: Could not find Fortran compiler"
	    echo "Plese export 'FC' variable before running this script"
	    exit 1
	fi
    fi
fi

# This is a very crude heuristic
comp=$(basename "$FC")
if [[ "$comp" = ifx || "$comp" = ifort || "$comp" = mpiifort ]]; then
    FCTYPE=intel
elif [[ "$comp" =~ ^gfortran || "$($FC --version 2>/dev/null | head -1)" =~ ^GNU ]]; then
    FCTYPE=gnu
else
    FCTYPE=unknown
fi

# Set default FFLAGS
if [[ $FCTYPE = intel ]]; then
    DEFAULT_FFLAGS="-i4 -check bounds,nouninit,noarg_temp_created -g -traceback -O0 -heap-arrays"
elif [[ $FCTYPE = gnu ]]; then
    warning_flags="-Wall -Wno-unused-dummy-argument -Wno-integer-division -Wno-unused-function -Wno-maybe-uninitialized"
    # TODO: We should probably trap also "invalid" floating-point exception,
    # but need to fix some tests first that currently trigger them :-(
    # We could also trap "underflow" but that seems a bit dangerous,
    # as it could crash simulations that would otherwise be fine?
    DEFAULT_FFLAGS="-Og -g --check=all -ffpe-trap=zero,overflow -fimplicit-none $warning_flags"
fi
# Existing FFLAGS are appended and thus should thake precedence
FFLAGS="${DEFAULT_FFLAGS} ${FFLAGS-}"

echo "Compiling with '$ESP' interface"
echo "Using $FC compiler"
$FC --version 2> /dev/null

cat << EOF > $VARS
ESP = $ESP
FC = $FC
LD = $FC
STDLIB_LIBDIR = $STDLIB_LIBDIR
STDLIB_MODULEDIR = $STDLIB_MODULEDIR
FCTYPE = $FCTYPE
FFLAGS = $FFLAGS
KERNEL = $(uname -s)
EOF
echo "Wrote this to '$VARS'"
echo "====================="
cat $VARS
