# we need the machine type to decide what to call things and what flags to use
MACHINE = $(shell uname -s)
# we need Matlab's path to know where to find includes and libraries
MATLABPATH = $(shell ./getmatlabpath.sh)

# here are the settings for Mac machines
# XXX ideally we would get this from 'mex', but that has so much bloat...
ifeq ($(MACHINE),Darwin)
  ENDING = mexmaci64
  LIBNAME = maci64
  LDFLAGS = -bundle
  # for the clang compilers we need to enable C++11
  ISCLANG=$(shell $(CXX) --version 2>/dev/null | grep -c clang)
  ifeq ($(ISCLANG),1)
    CXXFLAGS += -std=c++11
  endif
else ifeq ($(MACHINE),Linux)
  ENDING = mexa64
  LIBNAME = glnxa64
  CXXFLAGS += -fPIC
  LDFLAGS = -Wl,-rpath-link,$(MATLABPATH)/bin/glnxa64 -pthread -shared
else
    $(error Unrecognized machine type: $(MACHINE).)
endif

CXXFLAGS += -DNDEBUG -I$(MATLABPATH)/extern/include -O3
LDFLAGS += -L$(MATLABPATH)/bin/$(LIBNAME) -lmx -lmex -lmat -O3

MEX_BASE=calcseqw_cpp
MEX_NAME=$(MEX_BASE).$(ENDING)
OBJECTS=$(MEX_BASE).o calcobject.o

all: $(MEX_NAME)

$(MEX_NAME): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

calcobject.o: calcobject.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

${MEX_BASE}.o: ${MEX_BASE}.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	-rm -f $(MEX_NAME)

distclean: clean
	-rm -f $(OBJECTS)

.PHONY: all clean distclean
