#!/bin/sh

TMPFILE=/tmp/$$.out

FCTS_ROOT=${FCTS_ROOT:-./fcts}
MZN=${MZN:-./mistral-fz}
CANON=${CANON:-solns2dzn -c}

test_dir()
{
    DIR=$1
    for i in `find $FCTS_ROOT/$DIR -name "*.fzn"`; do
	rm -f $TMPFILE
	FDIR=`dirname $i`
	OPT=$FDIR/`basename $i .fzn`.opt
	EXP=$FDIR/`basename $i .fzn`.exp
	$MZN `cat $OPT` $i | $CANON 2>&1 > $TMPFILE
	if diff -uw $TMPFILE $EXP 2>&1 > /dev/null; then
	    echo $i pass
	else
	    echo $i failed
	    diff -uw $TMPFILE $EXP
	fi
    done
}

if test $# -eq 0; then
    test_dir builtins/int
else
    while test $# -gt 0
    do
	TESTDIR=$1
	shift
	test_dir $TESTDIR
    done
fi
