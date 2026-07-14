#
# Makefile to compile fms object files and assemble a library.
#
# This makefile does NOT work on its own -
# It needs to be invoked via the top-level Makefile
#
######################
# TARGET - what we're here for
####################
TARGET=libfms.a

##################
# Machine Dependencies
##################
MAKEFILE=Makefile

CONFIG=../CONFIGFMS

MODULEDIR=modules/

.SUFFIXES:

SUFFIXES=.f .F .f90 .F90

.PHONY: FORCE modules list clean

#################################
### Sources Headers and Objects #
#################################

all:		modules $(TARGET)

SRCS		= $(filter-out build_info.F90 openfms.F90,$(foreach suffix,$(SUFFIXES),$(wildcard *$(suffix))))
MOD_SRCS	= $(foreach suffix,$(SUFFIXES),$(wildcard $(MODULEDIR)/*$(suffix)))

INTFILES        = $(sort $(wildcard *.int))

OBJS		= $(addsuffix .o,$(basename $(SRCS)))
MOD_OBJS	= $(addsuffix .o,$(basename $(MOD_SRCS)))


$(TARGET): build_info.o openfms.o
	@echo "Creating FMS object library: $(TARGET)"
	@ar -crus $(TARGET) $(OBJS) $(MOD_OBJS)
ifeq ($(KERNEL),Darwin)
	ranlib -c $(TARGET)
endif

MODFILES = modules/globalmodule.mod modules/fmsmodule.mod modules/eispackmodule.mod \
	modules/bundlemodule.mod modules/bundleiomodule.mod modules/bundlecalcsmodule.mod \
	modules/spawnmodule.mod modules/elecstrucmodule.mod modules/randommodule.mod \
	modules/initialmodule.mod modules/samplingmodule.mod modules/particlemodule.mod \
	modules/trajectorymodule.mod modules/trajectoryiomodule.mod modules/trajectorycalcsmodule.mod \
	modules/restartmodule.mod modules/templatemodule.mod modules/toymodelmodule.mod\
	modules/qm_mm_module.mod modules/terachemmodule.mod modules/verletmodule.mod \
	modules/quanticsmodule.mod modules/rattlemodule.mod modules/thermomodule.mod \
	modules/smdmodule.mod modules/constraintmodule.mod modules/propagationmodule.mod \
	modules/minsearchmodule.mod modules/selectionmodule.mod modules/overlapmodule.mod

#################################
#  Module dependencies          #
#################################

include moddepend

moddepend: $(SRCS) ./Makemoddepend.pl
	@echo "Generating module dependencies"
	@./Makemoddepend.pl $(SRCS) openfms.F90

########################
###   Rules      #######
########################

FC_VERSION:=$(shell $(FC) --version 2>/dev/null | sed -n '1p')
BUILD_DATE:=$(shell date '+%Y-%m-%d %H:%M:%S %z')

build_info.F90: FORCE
	@OPENFMS_SOURCE_DIR=.. \
	OPENFMS_VERSION='v0.1-alpha' \
	OPENFMS_BACKEND='$(ESP)' \
	OPENFMS_BUILD_SYSTEM='Make' \
	OPENFMS_BUILD_TYPE='Makefile' \
	OPENFMS_FC='$(FC)' \
	OPENFMS_FC_ID='$(FC_VERSION)' \
	OPENFMS_FC_VERSION='' \
	OPENFMS_FFLAGS='$(FFLAGS)' \
	../generate_build_info.sh $@

# Always recompile build_info.F90 to keep the git revision up to date.
build_info.o: build_info.F90 openfms.o $(OBJS) $(MOD_OBJS)
	$(FC) -c $(FFLAGS) -I $(MODULEDIR) -DOPENFMS_BUILD_DATE="'$(BUILD_DATE)'" $<

%.o : %.F
	@echo "Compiling $<"
	@$(FC) -c $(FFLAGS) -I $(MODULEDIR) $<

%.o : %.f
	@echo "Compiling $<"
	@$(FC) -c $(FFLAGS) -I $(MODULEDIR) $<

%.o : %.f90
	@echo "Compiling $<"
	@$(FC) -c $(FFLAGS) -I $(MODULEDIR) $<

%.o : %.F90
	@echo "Compiling $<"
	@$(FC) -c $(FFLAGS) -I $(MODULEDIR) $<

$(MODFILES): modules
# NOTE: It's VERY important to have the recipe for MOD_OBJS empty (note the semicolon at the end),
# otherwise the modules/*.o files will get implicit rules and will be compiled from this Makefile,
# while they should be compiled by modules/Makefile
# See: https://www.gnu.org/software/make/manual/html_node/Empty-Recipes.html
$(MOD_OBJS): modules ;

modules: $(MAKEFILE) $(CONFIG)
	@cd $(MODULEDIR); $(MAKE) -r

# Rebuild everything if modules have changes
$(OBJS): $(MAKEFILE) $(CONFIG)

########################
### Additional Targets #
########################
list:
	@echo "Sources: " $(SRCS)
	@echo "Objects: " $(OBJS)
	@echo "Module sources: " $(MOD_SRCS)

clean:
	$(RM) $(TARGET) moddepend *.o *.mod build_info.F90
	@cd $(MODULEDIR); make clean
