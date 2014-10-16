Installing Numberjack
=====================

Numberjack offers the ability to use several high-performance solvers, some of
which are required to be installed separately or have their own dependencies.
The source-code for building Mistral, `Toulbar2
<https://mulcyber.toulouse.inra.fr/projects/toulbar2/>`_, `Minisat
<http://minisat.se>`_, and `Walksat
<http://www.cs.rochester.edu/u/kautz/walksat/>`_ is included in the Numberjack
distribution and interfaces to additional solvers are available.

To install Numberjack from source simply run `python setup.py build` from the
Numberjack directory, followed by `python setup.py install`. This requires the
following to be installed:

* python-dev
* swig
* libxml2-dev
* zlib1g-dev 
* libgmp-dev


Building additional solver interfaces
-------------------------------------

The source-code for building Mistral, Toulbar2, Minisat, and Walksat is included
with the Numberjack distribution and interfaces to additional solvers are
available.

The following solvers are entirely optional and are not required to use
Numberjack however they can be a valuable tool to have. If any of the solvers
cannot be found on the system, then their interface will be disabled within
Numberjack. If you add one of these solvers after already installing Numberjack,
please reinstall Numberjack to enable the new interface.


Gurobi Optimizer
^^^^^^^^^^^^^^^^

Numberjack provides an interface to the mathematical programming solver, Gurobi.
To use it in Numberjack, `Download and install Gurobi Optimizer
<http://www.gurobi.com/download/gurobi-optimizer>`_ from their website.
Numberjack has been tested against Gurobi versions 5.1, 5.5, and 5.6.

Numberjack will try to automatically find the Gurobi install location by first
checking the :envvar:`GUROBI_HOME` environment variable which should be set
according to the Gurobi install guide. If not, it will try to find the path
based on the location of the `gurobi_cl` executable.

.. code-block:: bash

    export GUROBI_HOME="/opt/gurobi550/linux64"    # Example Linux install dir
    export GUROBI_HOME="/Library/gurobi550/mac64"  # Example Mac OS X install dir



CPLEX
^^^^^

Numberjack provides an interface to `IBM ILOG CPLEX Optimizer
<http://www.ibm.com/software/commerce/optimization/cplex-optimizer/>`_ and has
been tested against CPLEX versions 12.5, and 12.5.1.

Numberjack will try to automatically find the CPLEX install location by first
checking the environment variable :envvar:`CPLEXDIR`. If this is not set it will try to
find it based on the location of the `cplex` executable.

.. code-block:: bash

    export CPLEXDIR="/opt/ibm/ILOG/CPLEX_Studio1251/cplex"
    export CPLEXDIR="/Applications/IBM/ILOG/CPLEX_Studio1251/cplex"



SCIP
^^^^

SCIP is an open-source MIP solver available from: http://scip.zib.de, please
download version 3.1.0. The simplest way to install SCIP is to `download the
source code of SCIP Optimization Suite v3.1.0
<http://scip.zib.de/download.php?fname=scipoptsuite-3.1.0.tgz>`_ and place
scipoptsuite-3.1.0.tgz in the same folder as Numberjack, alongside this file, so
that it can be compiled with the necessary flags.

If you would like to compile SCIP yourself from a different location, then
please set the environment variable :envvar:`ZIBPATH` to the scipoptsuite
directory and compile the static library using the following options:

.. code-block:: bash

    export ZIBPATH=path_to/scipoptsuite-3.1.0
    make scipoptlib ZIMPL=false ZLIB=false READLINE=false GAMS=false GMP=false LEGACY=true SPX_LEGACY=true
