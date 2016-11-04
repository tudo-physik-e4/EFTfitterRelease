# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := -lMinuit $(shell root-config --libs)

# compiler and flags
CXX          = g++
CXXFLAGS     =  -g -Wall -fPIC -Wno-deprecated -O2
LD           = @LD@
LDFLAGS      =  -g -O2
SOFLAGS      = -shared

# standard commands
RM           = rm -f
MV           = mv
ECHO         = echo

# add ROOT flags
CXXFLAGS    += $(ROOTCFLAGS) 

# ----------------------------------------------------------------------

# The following definitions depend on the setup of the system where
# the project is being compiled. If BAT is installed in the standard
# system search path or the installation directory is defined in the
# BATINSTALLDIR environmental variable then the lines below are correct
# and the compilation will work

CXXFLAGS    += $(shell bat-config --cflags)
LIBS        += $(shell bat-config --libs)

# List of all classes (models) used in the program
CXXSRCS      = \
	AnomCouplings.cxx \
	model_specific/XsecCoeffs.cxx \
	model_specific/AccClass.cxx \
	../../src/tinyxml2/tinyxml2.cxx \
	../../src/mvclocal/BCMVCPhysicsModel_EFT.cxx \
	../../src/mvclocal/BCMVCDataModel.cxx \
	../../src/mvclocal/BCMVCMeasurement.cxx \
	../../src/mvclocal/BCMVCObservable.cxx \
	../../src/mvclocal/BCMVCombination_EFT.cxx \
	../../src/mvclocal/BCMVCUncertainty.cxx \

# ----------------------------------------------------------------------

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXXSRCS))
EXEOBJS      =
MYPROGS     = \
	runAnomCouplings 

GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d $(MYPROGS)


# targets
all : directories myfit

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

directories:
	mkdir -p results

myfit : runAnomCouplings.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) runAnomCouplings.o $(CXXOBJS) $(LIBS) -o runAnomCouplings

print :
	echo compiler  : $(CXX)
	echo c++ srcs  : $(CXXSRCS)
	echo c++ objs  : $(CXXOBJS)
	echo c++ flags : $(CXXFLAGS)
	echo libs      : $(LIBS)
	echo so flags  : $(SOFLAGS)

	echo rootlibs  : $(ROOTLIBS)

