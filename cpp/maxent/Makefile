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

MEX_EVAL_BASE=cppevalmaxent
MEX_EVAL_NAME=$(MEX_EVAL_BASE).$(ENDING)

MEX_SIM_BASE=cppsimmaxent
MEX_SIM_NAME=$(MEX_SIM_BASE).$(ENDING)

MEX_STATS_BASE=cppslowstats
MEX_STATS_NAME=$(MEX_STATS_BASE).$(ENDING)

OBJECTS_EVAL=$(MEX_EVAL_BASE).o couplings.o evalobject.o evaluator.o utils.o
OBJECTS_SIM=$(MEX_SIM_BASE).o couplings.o cppsimmaxent.o evaluator.o matlab_random.o mcmc.o simobject.o utils.o
OBJECTS_STATS=$(MEX_STATS_BASE).o couplings.o cppslowstats.o evaluator.o slowstatsobject.o utils.o

all: evaluator simulator slowstats

evaluator: $(MEX_EVAL_NAME)

simulator: $(MEX_SIM_NAME)

slowstats: $(MEX_STATS_NAME)

$(MEX_EVAL_NAME): $(OBJECTS_EVAL)
	$(CXX) $(LDFLAGS) -o $@ $^

# the simulation uses the undocumented utIsInterruptPending function
# we thus need to link to libut
$(MEX_SIM_NAME): $(OBJECTS_SIM)
	$(CXX) $(LDFLAGS) -lut -o $@ $^

$(MEX_STATS_NAME): $(OBJECTS_STATS)
	$(CXX) $(LDFLAGS) -o $@ $^

couplings.o: couplings.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

evalobject.o: evalobject.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

evaluator.o: evaluator.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

utils.o: utils.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

matlab_random.o: matlab_random.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mcmc.o: mcmc.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

simobject.o: simobject.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

slowstatsobject.o: slowstatsobject.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

${MEX_EVAL_BASE}.o: ${MEX_EVAL_BASE}.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

${MEX_SIM_BASE}.o: ${MEX_SIM_BASE}.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

${MEX_STATS_BASE}.o: ${MEX_STATS_BASE}.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	-rm -f $(MEX_EVAL_NAME) $(MEX_SIM_NAME) $(MEX_STATS_NAME)

distclean: clean
	-rm -f $(OBJECTS_EVAL) $(OBJECTS_SIM) $(OBJECTS_STATS)

.PHONY: all evaluator simulator slowstats clean distclean
