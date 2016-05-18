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

CXXFLAGS    += -I. -I./include -I$(BATINSTALLDIR)/include
LIBS        += -L$(BATINSTALLDIR)/lib -lBATmodels -lBATmvc -lBAT $(ROOTLIBS)

# List of all classes (models) used in the program
CXXSRCS      = \
	EFTfitter.cxx \
	../../src/tinyxml2/tinyxml2.cxx \
	../../src/mvclocal/BCMVCPhysicsModel.cxx \
	../../src/mvclocal/BCMVCDataModel.cxx \
	../../src/mvclocal/BCMVCMeasurement.cxx \
	../../src/mvclocal/BCMVCObservable.cxx \
	../../src/mvclocal/BCMVCombination_EFT.cxx \
	../../src/mvclocal/BCMVCUncertainty.cxx \


# ----------------------------------------------------------------------

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXXSRCS))
EXEOBJS      =
MYPROGS     = \
	runEFTfitter

GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d $(MYPROGS)


# targets
all : myfit

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

myfit : runEFTfitter.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) runEFTfitter.o $(CXXOBJS) $(LIBS) -o runEFTfitter

print :
	echo compiler  : $(CXX)
	echo c++ srcs  : $(CXXSRCS)
	echo c++ objs  : $(CXXOBJS)
	echo c++ flags : $(CXXFLAGS)
	echo libs      : $(LIBS)
	echo so flags  : $(SOFLAGS)

	echo rootlibs  : $(ROOTLIBS)

