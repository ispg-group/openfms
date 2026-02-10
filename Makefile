#############
### Targets #
#############

#Configuration information read from CONFIGFMS:
include CONFIGFMS
export FC
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

.PHONY: makefmslib test testclean updatetestref clean veryclean list

MODULEDIR=Modules/

SUFFIXES=.f .c .F

####################
# How to build libfms.a
####################

FMSLIB = src/libfms.a

$(FMSLIB): makefmslib

# NOTE: It's important to use $(MAKE) for working parallel compilation!
# https://stackoverflow.com/a/60706372/3682277
makefmslib: CONFIGFMS
	@cd src; $(MAKE) -r

####################
# Linking the executable
####################
bin:
	mkdir bin

# This is standalone OpenFMS, without external electronic structure
bin/$(PROGBASE).zero: $(FMSLIB) $(MAKEFILE) bin
		@rm -f bin/$(PROGBASE)
		@echo;echo "Linking $(PROGRAM) ..."
		$(LD) $(FMSLIB) -I$(MODULEDIR) -o $(TARGET) $(LIBS) $(LDFLAGS)
		@ln -s $(PROGRAM) bin/$(PROGBASE)
		@echo "done"

# OpenFMS with TeraChem interface
bin/$(PROGBASE).tc:  $(FMSLIB) $(MAKEFILE) bin
		@rm -f bin/$(PROGBASE)
		@echo;echo "Linking $(PROGRAM) ..."
		$(LD) $(FFLAGS) $(FMSLIB) -o $(TARGET) $(LIBS) $(LDFLAGS)
		@ln -s $(PROGRAM) bin/$(PROGBASE)
		@echo "done"

bin/$(PROGBASE).quantics: $(FMSLIB) $(MAKEFILE) bin
		@rm -f bin/$(PROGBASE)
		@echo;echo "Linking $(PROGRAM) ..."
		$(LD) $(FMSLIB) -I$(MODULEDIR) -o $(TARGET) $(LIBS) $(LDFLAGS) -J$(QUANTICS_OBJ)/include -L$(QUANTICS_DIR)/bin/dyn_libs -lsrf -lusrf -lsqlite3 -Wl,-rpath=$(QUANTICS_DIR)/bin/dyn_libs
		@ln -s $(PROGRAM) bin/$(PROGBASE)
		@echo "done"

# End-To-End tests
test: ${TARGET}
	tests/test.sh "${TARGET}" "$(TEST)" test

# Unit tests
unittest: makefmslib CONFIGFMS
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

veryclean: clean
	rm bin/*
	rm CONFIGFMS
