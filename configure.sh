#!/bin/bash
# Embarassigly barebones configure.sh, which generates the CONFIGFMS file used by make
# The first parameter to this script determines the interface that is compiled.
# Currently, only the MPI-based interface to TeraChem is available.
# By default, we only compile with the model potentials ('zero').

# The second parameter to this script specifies the location of the stdlib libraries (optional).
# Unfortunately, this means that if you need to specify the location of the stdlib libraries, you also need to 
# specify the type of electronic-structure interface. Use 'zero' for this if you do not want any interface 
# and only models potentials.

# Output file name
VARS=CONFIGFMS

# Electronic structure interface
ESP=${1-zero}
if [[ $ESP != "tc" && $ESP != "zero" && $ESP != "quantics" ]]; then
    echo "ERROR: Unsupported interface target '$ESP'"
    echo "Supported targets: 'zero', 'tc'"
    exit 1
fi

# Fortran's STDLIB package
STDLIB_LIBDIR=${2-none}
STDLIB_MODULEDIR=${3-none}
if [[ $STDLIB_LIBDIR != "none" && $STDLIB_MODULEDIR = "none" ]]; then
	echo "ERROR: Please provide location of the STDLIB module directory"
	exit 1
fi
if [[ $STDLIB_LIBDIR = "none" ]]; then
	STDLIB_LIBDIR=""
	STDLIB_MODULEDIR=""
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
