[![codecov](https://codecov.io/gh/ispg-group/openfms/graph/badge.svg?token=H1IEDF52ZX)](https://codecov.io/gh/ispg-group/openfms)

# OpenFMS: ab initio multiple spawning code

This is an open-source version of the FMS90 code by Todd Martínez et al

The code is based on the `FMS-TeraChem-Simple` branch, and has been significantly pruned and refactored.
It currently contains the MPI interface to TeraChem and various analytical (toy) potentials.

## Installation

To compile the code and run the test suite, the following sequence should work on a typical unix machine with GCC toolchain intalled.

```console
./configure.sh
make -j
make unittest
make test
```

Note: The `configure.sh` script is very barebones, you might need to tweak the CONFIGFMS file that it generates.

In the default (`zero`) compilation target, only the analytical model potentials are available.

To compile with TeraChem interface, run:

```console
./configure.sh tc
```

In this case, the compiler needs to be changed to `mpif90` using the MPI implementation
that is matching the version that was used to compile TeraChem (typically MPICH or IntelMPI).
To compile MPICH, we recommend using the included `.github/install_mpich.sh` script.
Note that MPICH version >=4.0 is strongly encouraged as earlier versions
contained bugs that were affecting the TeraChem interface.

## Development

Please see the [Contribution Guide](CONTRIBUTING.md).
