#!/bin/sh

LOCAL_LIB=`pwd`/local_lib
FZN_LIB=`pwd`/fzn

echo  '  Remember to add the following lines to your .profile or .bashrc:'
echo  '    export PYTHONPATH='$LOCAL_LIB':$PYTHONPATH'
echo  '    export PATH='$LOCAL_LIB':$PATH:'$FZN_LIB
echo  '  then run "source ~/.profile" or "source ~/.bashrc"'
