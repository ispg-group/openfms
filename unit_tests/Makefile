# All source files in this folder should have .F90 suffix!
SUFFIXES = .F90
LIBFMS = ../src/libfms.a
INCFLAGS = -I../src/modules
MAKEFILE = Makefile
CONFIG = ../CONFIGFMS

.PHONY: clean run

# Pattern rules
%.o : %.F90
	@echo "Compiling $<"
	@$(FC) -Wall $(FFLAGS) $(INCFLAGS) -c $<

TESTFILES=$(wildcard test_*.F90) $(wildcard test_*.f90)
TESTOBJECTS=$(addsuffix .o,$(basename $(TESTFILES)))
TESTDRIVEOBJS = testdrive.o testutils.o

# Default target, run the tests
run: test
	./test

# Create a driver binary "test"
test: main.o $(TESTOBJECTS) $(TESTDRIVEOBJS)
	$(FC) $(LDFLAGS) main.o $(TESTDRIVEOBJS) $(TESTOBJECTS) $(LIBFMS) -o test

main.o: main.F90 $(TESTOBJECTS) $(TESTDRIVEOBJS)

$(TESTOBJECTS): $(LIBFMS) $(TESTDRIVEOBJS)

testutils.o: testutils.F90 testdrive.o

# Rebuild everything if Makefile or CONFIGFMS changes
testdrive.o: $(MAKEFILE) $(CONFIG)

clean:
	rm -f *.mod *.o ./test

