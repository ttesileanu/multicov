#!/bin/bash

CXX="g++"

# figure out the paths
# XXX this is very non-portable! Should run mex instead!
echo "Finding matlab mex path..."
MACHINE=`uname -s`
case $MACHINE in
Darwin)
  MATLABPATH=`echo /Applications/MATLAB_R20*.app`
  ENDING="mexmaci64"
  LIBNAME="maci64"
  CXXFLAGS=""
  LDOPTIONS="-bundle"
  CXX="clang++"
  command -v ${CXX} >/dev/null 2>&1 || CXX="g++"
  (${CXX} --version 2>/dev/null | grep clang) > /dev/null; tmp=$?
  if [ "$tmp" -eq 0 ] ; then
    CXXFLAGS="-std=c++11 $CXXFLAGS"
  fi
  ;;
Linux)
  # this is more complicated than it needs to be because matlab is silly and
  # adds some weird character at the beginning and end of its output...
  if [ ! -e "`which matlab`" ] ; then
    echo "Couldn't find matlab executable."
    exit 1
  fi
  MATLABPATH=`matlab -nosplash -nodesktop -r "disp(matlabroot); exit" | tail -n2 | head -n1`
  ENDING="mexa64"
  LIBNAME="glnxa64"
  CXXFLAGS="-fPIC"
  LDOPTIONS="-Wl,-rpath-link,$MATLABPATH/bin/glnxa64 -pthread -shared"
  ;;
*)
  echo "Unrecognized machine."
  exit 1
  ;;
esac

if [ ! -d "$MATLABPATH" ] ; then
  echo "Couldn't find matlab root path."
  exit 1
fi

echo "Compiling calcseqw_cpp..."
rm -f calcobject.o calcseqw_cpp.o
${CXX} -c -DNDEBUG ${CXXFLAGS} -I${MATLABPATH}/extern/include -O3 -o calcobject.o calcobject.cc

${CXX} -c -DNDEBUG ${CXXFLAGS} -I${MATLABPATH}/extern/include -O3 -o calcseqw_cpp.o calcseqw_cpp.cc

echo "Linking calcseqw_cpp..."
${CXX}  ${LDOPTIONS} -L${MATLABPATH}/bin/${LIBNAME} -lmx -lmex -lmat -O3 -o calcseqw_cpp.${ENDING} calcseqw_cpp.o calcobject.o
