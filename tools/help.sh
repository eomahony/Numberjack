#!/bin/sh

LOCAL_LIB=`pwd`/local_lib

echo  '  Remember to add the following lines to your .profile or .bashrc:\n'
echo  '    PYTHONPATH='$LOCAL_LIB':$PYTHONPATH'
echo  '    export PYTHONPATH\n'
echo  '  and'
echo  '    PATH='$LOCAL_LIB':$PATH'
echo  '    export PATH\n'
echo  '  then run "source ~/.profile" or "source ~/.bashrc"'
