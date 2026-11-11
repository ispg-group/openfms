#############
### Targets #
#############

#Configuration information read from CONFIGFMS:
include CONFIGFMS
export FC
export FCTYPE
export LD
export FFLAGS
export LDFLAGS
export KERNEL
export ESP
export MV = mv
export SHELL = bash
export CP = cp
export RM = rm -f

# By default run all end-to-end tests in tests/
TEST=all

PROGBASE	= openfms
PROGRAM		= $(PROGBASE).$(ESP)
TARGET		= bin/$(PROGRAM)
MAKEFILE	= Makefile
QML_path := $(CURDIR)/QuantumModelLib

# Decide what to compile
ifeq ($(ESP),zero)
   all:	$(TARGET)
endif
ifeq ($(ESP),tc)
   all: $(TARGET)
   FFLAGS += -DTeraChem
endif
ifeq ($(ESP),quantics)
   all:	$(TARGET)
   FFLAGS += -DQuantics
   QUANTICS_OBJ = $(QUANTICS_DIR)/object/x86_64/gfortran/
   LIBS = $(QUANTICS_OBJ)/versions.o
   LIBS += ${QUANTICS_OBJ}/quanticslib.a
   LIBS += ${QUANTICS_OBJ}/opfuncs.a
   LIBS += ${QUANTICS_OBJ}/quanticsmod.a
   LIBS += ${QUANTICS_OBJ}/includes.a
   LIBS += ${QUANTICS_OBJ}/globinc.a
   LIBS += ${QUANTICS_OBJ}/libnum.a
   LIBS += ${QUANTICS_OBJ}/libsys.a
   LIBS += ${QUANTICS_OBJ}/liblapack.a
   LIBS += ${QUANTICS_OBJ}/libblas.a
   LIBS += ${QUANTICS_OBJ}/libomp.a
endif
ifeq ($(ESP),qml)
   all: $(TARGET)
   FFLAGS += -DQuantumModelLib
   LIBS += $(QML_path)/libQMLFull.a
endif
.PHONY: FORCE test testclean updatetestref clean veryclean list

MODULEDIR=Modules/

SUFFIXES=.f .c .F

####################
# How to build libfms.a
####################

FMSLIB = src/libfms.a
MAIN = src/openfms.o src/build_info.o

$(MAIN): $(FMSLIB)
# NOTE: It's important to use $(MAKE) for working parallel compilation!
# https://stackoverflow.com/a/60706372/3682277
$(FMSLIB): CONFIGFMS FORCE bin
ifeq ($(ESP),qml)
	@cd $(QML_path); $(MAKE) -r LAPACK=0 OPT=1 OMP=1 REAL=64 INT=4
endif
	@cd src; $(MAKE) -r
####################
# Linking the executable
####################
bin:
	mkdir bin

# This is standalone OpenFMS, without external electronic structure
bin/$(PROGBASE).zero: $(FMSLIB) $(MAIN) $(MAKEFILE)
		@echo;echo "Linking $(PROGRAM) ..."
		$(LD) $(MAIN) $(FMSLIB) -I$(MODULEDIR) -o $(TARGET) $(LIBS) $(LDFLAGS)
		@rm -f bin/$(PROGBASE)
		@ln -s $(PROGRAM) bin/$(PROGBASE)
		@echo "done"

# OpenFMS with TeraChem interface
bin/$(PROGBASE).tc: $(FMSLIB) $(MAIN) $(MAKEFILE)
		@rm -f bin/$(PROGBASE)
		@echo;echo "Linking $(PROGRAM) ..."
		$(LD) $(FFLAGS) $(MAIN) $(FMSLIB) -o $(TARGET) $(LIBS) $(LDFLAGS)
		@ln -s $(PROGRAM) bin/$(PROGBASE)
		@echo "done"

bin/$(PROGBASE).quantics: $(FMSLIB) $(MAIN) $(MAKEFILE)
		@rm -f bin/$(PROGBASE)
		@echo;echo "Linking $(PROGRAM) ..."
		$(LD) $(MAIN) $(FMSLIB) -I$(MODULEDIR) -o $(TARGET) $(LIBS) $(LDFLAGS) -J$(QUANTICS_OBJ)/include -L$(QUANTICS_DIR)/bin/dyn_libs -lsrf -lusrf -lsqlite3 -Wl,-rpath=$(QUANTICS_DIR)/bin/dyn_libs
		@ln -s $(PROGRAM) bin/$(PROGBASE)
		@echo "done"

# OpenFMS with QuantumModelLib interface
bin/$(PROGBASE).qml: $(FMSLIB) $(MAIN) $(MAKEFILE)
		@rm -f bin/$(PROGBASE)
		@echo;echo "Linking $(PROGRAM) ..."
		$(LD) $(MAIN) $(FMSLIB) -o $(TARGET) $(LIBS) -fopenmp $(LDFLAGS)
		@ln -s $(PROGRAM) bin/$(PROGBASE)
		@echo "done"

# End-To-End tests
test: ${TARGET}
	tests/test.sh "${TARGET}" "$(TEST)" test

# Unit tests
unittest: $(TARGET)
	@cd unit_tests; $(MAKE) -r

# Clean all test folders.
testclean:
	tests/test.sh "${TARGET}" "$(TEST)" clean

# This will automatically generate new reference data for a given test
updatetestref: ${TARGET}
	tests/test.sh "${TARGET}" "$(TEST)" makeref

########################
### Additional Targets #
########################
list:
	@cd src; make list

clean:
	@cd src; make clean
ifeq ($(ESP),qml)
	@cd $(QML_path); make clean
endif
veryclean: clean
	rm bin/*
	rm CONFIGFMS
