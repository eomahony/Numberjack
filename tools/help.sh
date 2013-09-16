#!/bin/sh

LOCAL_LIB=`pwd`/local_lib

echo  '  Remember to add the following lines to your .profile or .bashrc:'
echo  '    export PYTHONPATH='$LOCAL_LIB':$PYTHONPATH'
echo  '    export PATH='$LOCAL_LIB':$PATH'
echo  '  then run "source ~/.profile" or "source ~/.bashrc"'
