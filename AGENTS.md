# AGENTS.md

OpenFMS is an open-source Fortran implementation of ab initio multiple spawning (AIMS / FMS90, Martínez et al.), with an MPI interface to TeraChem and a set of analytical (toy) model potentials.

## Build and test

```console
./configure.sh          # writes CONFIGFMS (default 'zero' backend = toy models only)
./configure.sh tc       # TeraChem backend; needs mpifort/mpiifort built with the same MPI as TeraChem
make -j                 # builds src/libfms.a, links bin/openfms.<backend> (+ symlink bin/openfms)
make unittest           # builds and runs unit_tests/test
make test               # all end-to-end tests in tests/
make test TEST=AIMS     # a single end-to-end test (directory name under tests/)
make updatetestref TEST=AIMS   # overwrite existing *.ref files for ONE test with current output
make testclean
```

- `CONFIGFMS` is generated; tweak `FC`/`FFLAGS` there or export `FC`/`FFLAGS` before `./configure.sh` (exported `FFLAGS` are appended to the debug defaults: `-Og -g --check=all -ffpe-trap=zero,overflow -fimplicit-none`).
- Run a single unit-test suite: `unit_tests/test <suite>` (suite names: `ParticleModule`, `TrajectoryModule`, `BundleModule`, `OverlapModule`, as registered in `unit_tests/main.F90`).
- CMake (zero backend only): `cmake -S . -B build -DBUILD_TESTING=ON && cmake --build build -j && ctest --test-dir build --output-on-failure`. Source lists live in `cmake/*.cmake` and must be kept in sync with the Makefile build when adding files.
- Module dependencies for Make are generated into `src/moddepend` by `src/Makemoddepend.pl`.

## Lint and format

Hooks are run via `prek` (config in `prek.toml`): `prek run -a` (all hooks) or `prek run -a format` (formatter only). Hooks include:
- `fortitude` linter (`fortitude.toml`, line length 132, single-quoted strings),
- `fprettify` formatter via `./autoformat.py` (`.fprettify.rc`: 3-space indent, C-style relational operators `< > == /=`),
- `./check_module_usage.sh`,
- `ruff` for the Python helper scripts.

Developer docs (FORD): `ford ford_config.md` → `developer_docs/`.

## Code conventions

- `implicit none` everywhere; every dummy argument has an `intent`.
- Use `real(DefReal)`, `integer(DefInt)`, `complex(DefComp)` kinds from `GlobalModule`.
- Modules are `private` by default with explicit `public` lists.
- `use SomeModule, only: ...` is enforced, except the modules in `allow-bare-use` in `fortitude.toml`.
- **Never** use `only:` with `ParticleModule`, `TrajectoryModule` or `BundleModule`: they overload `assignment(=)` for `t_Particle`, `t_Trajectory`, `t_TrajectoryBundle`, and `only:` lists silently drop the overload (enforced by `check_module_usage.sh`).
- Errors are raised via `FMS_DieError` (error handler is swappable via `set_error_handler`, which the unit tests use).

## Architecture

- `src/openfms.F90` is the main program: reads `Control.dat` (`FMS_ReadNameList` in `src/read_namelist.f90`, a single `&control` namelist), seeds the PRNG, initialises electronic structure (`FMS_ESInit`) and the interface, builds the trajectory bundle, sets initial conditions or restarts from `Last_Bundle.txt`, then loops `FMS_Dynamics` / `FMS_Output` until the simulation time is reached, and finally `FMS_Shutdown`.
- Core data model: `t_Particle` (atoms) → `t_Trajectory` (one Gaussian basis function / TBF on a given state) → `t_TrajectoryBundle` (the full wavefunction). `*CalcsModule` and `*IOModule` hold computations and I/O for trajectories/bundles; `SpawnModule` handles spawning of new TBFs; `OverlapModule`, `PropagationModule`, `VerletModule` cover overlaps and integration.
- Electronic structure dispatch is in `ElecStrucModule` (`select case (gliModel)`): `FMSZERO`/toy models in `ToyModelModule` (selected by `iMethod`), `TEMPLATE` (`TemplateModule`, external programs via input/output templates), `TC` (`TerachemModule.F90`, compiled only with `-DTeraChem`), `QUANTICS` (`-DQuantics`). Backend-specific files use `.F90` for preprocessing.
- New namelist input variables must be declared and added to the `&control` namelist in `src/read_namelist.f90`.

## Tests

- End-to-end tests: each `tests/<UPPERCASE_DIR>/` contains `Control.dat` (and `Geometry.dat`, etc.) plus `*.ref` reference files; only directories starting with a capital letter are run (e.g. `tests/swarm` is skipped). `tests/test.sh` runs the binary (or a per-directory `test.sh` if present, e.g. `CMDLINE`, which writes `errors.dat`), then diffs each output against its `.ref`, tolerating small numerical noise via `tests/numdiff.py`. `errors.dat` must match exactly. Restart tests copy `Last_Bundle.txt.initial` first. New tests need their `.ref` files created manually; `updatetestref` only refreshes existing ones.
- Unit tests use the bundled `test-drive` framework (`unit_tests/testdrive.F90`). New suites must be named `unit_tests/test_<suite>.F90` (Makefile globs on this) and registered in both `unit_tests/main.F90` and `cmake/OpenFMSUnitTests.cmake`.
- CI (`.github/workflows/build.yml`) builds with GCC (incl. GCC 9 and macOS), Intel, Intel classic, AOCC, a sanitizer build (ASan/UBSan), CMake, and the TeraChem backend with MPICH/IntelMPI — keep code portable across these compilers.
