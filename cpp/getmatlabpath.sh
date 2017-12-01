#! /bin/bash

MACHINE=`uname -s`
case $MACHINE in
Darwin)
  MATLABPATH=`echo /Applications/MATLAB_R20*.app`
  ;;
Linux)
  if [ ! -e "`which matlab`" ] ; then
    echo "Couldn't find Matlab executable."
    exit 1
  fi
  MATLABPATH=`matlab -nosplash -nodesktop -r "disp(matlabroot); exit" | tail -n2 | head -n1`
  ;;
*)
  echo "Unrecognized machine."
  exit 1
  ;;
esac

if [ ! -d "$MATLABPATH" ] ; then
  echo "Couldn't find Matlab root path."
  exit 1
fi

echo $MATLABPATH
