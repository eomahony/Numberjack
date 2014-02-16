# Numberjack
## A Python constraint programming platform 

## What is Numberjack?

Numberjack is a modelling package written in Python for constraint programming.
Python benefits from a large and active programming community, Numberjack is
therefore a perfect tool to embed CP technology into larger applications. It is
designed to support a number of underlying C/C++ solvers seamlessly and
efficiently. Currently, there are six available back-ends: three mixed integer
programming solvers ([Gurobi][gurobiopt], [CPLEX][cplex], and [SCIP][scipopt]),
two satisfiability solvers ([MiniSat][minisat] and [Walksat][walksat]), and a
constraint programming solver (Mistral).

* Numberjack offers a high level constraint programming language
* Numberjack directly benefits from python's features and modules
* Numberjack uses efficient underlying C/C++ solvers.

## Installation

Numberjack offers the ability to use several high-performance solvers, some of
which are required to be installed separately or have their own dependencies.
The source-code for building Mistral, [Minisat][minisat], and [Walksat][walksat]
is included in the Numberjack distribution and interfaces to additional solvers
are available.

To install Numberjack from source simply run `make` from the Numberjack
directory, followed by `make install`. This requires the following to be
installed:

* python 2
* swig
* libxml2-dev
* zlib1g-dev
* python-dev
* libgmp-dev

[minisat]: http://minisat.se
[walksat]: http://www.cs.rochester.edu/u/kautz/walksat/


## Building the Solvers

### Gurobi Optimizer

Numberjack provides an interface to the mathematical programming solver,
[Gurobi][gurobiopt]. To use it in Numberjack, please follow the following steps:

1. [Download and install Gurobi Optimizer][gurobiopt] from their website.
   Numberjack has been tested against Gurobi versions 5.1, 5.5, and 5.6.

2. Numberjack will try to find the Gurobi, install location by first checking
   the `GUROBI_HOME` environment variable which should be set according to the
   Gurobi install guide. If not, it will try to find the path based on the
   location of the `gurobi_cl` executable.

   ```bash
   export GUROBI_HOME="/opt/gurobi550/linux64"    # Example Linux install dir
   export GUROBI_HOME="/Library/gurobi550/mac64"  # Example Mac OSX install dir
   ```

3. In the Numberjack directory, move the folder `available_interfaces/gurobi` to
   `solvers/gurobi`.

4. Run `make`, followed by `make install` from the Numberjack directory.

[gurobiopt]: http://www.gurobi.com/download/gurobi-optimizer
[gurobiqs]: http://www.gurobi.com/documentation/current/quick-start-guide/


### CPLEX

Numberjack provides an interface to [IBM ILOG CPLEX Optimizer][cplex]. If you
have it installed on your system, follow these steps to use in in Numberjack:

1. If the executable `cplex` is in the PATH Numberjack will be able to find the
   location of CPLEX on your system. If not, you should set the environment
   variable `CPLEXDIR` to the location where CPLEX is installed, for example:

    ```bash
    export CPLEXDIR="/opt/ibm/ILOG/CPLEX_Studio1251/cplex"
    export CPLEXDIR="/Applications/IBM/ILOG/CPLEX_Studio1251/cplex"
    ```

2. In the Numberjack directory, move the folder `available_interfaces/cplex` to
   `solvers/cplex`.

3. Run `make`, followed by `make install` from the Numberjack directory.

[cplex]: http://www.ibm.com/software/commerce/optimization/cplex-optimizer/


### SCIP
Due to licensing restrictions Numberjack cannot include the sources for SCIP.
The SCIP source code is available from:
[http://scip.zib.de/download.shtml][scipopt] (available for free under an
academic license). Please download version 3.0.1.

[scipopt]: http://scip.zib.de/download.shtml

1. Download the source code of SCIP Optimization Suite v3.0.1, unarchive it,
   however, do not install nor compile it.

    ```bash
    tar zxf scipoptsuite-3.0.1.tgz
    ```

2. In the scipoptsuite directory, unarchive the SCIP and SoPlex archives.

    ```bash
    cd scipoptsuite-3.0.1/
    tar zxf scip-3.0.1.tgz
    tar zxf soplex-1.7.1.tgz
    ```

3. Set an environment variable ZIBPATH to the scipoptsuite directory, for
   example:

    ```bash
    export ZIBPATH=path_to/scipoptsuite-3.0.1
    ```

4. Move the folder `available_interfaces/scip` to `solvers/scip`.

5. Run `make`, followed by `make install` from the Numberjack directory.


### Osi Solvers
Some of the solvers cannot be included in the source tree.
Before you get started with the compiling Osi solvers for Numberjack you must first get the Osi module working.

Download Osi from: [http://www.coin-or.org/download/source/Osi/][osi]

To do this you need to download Osi-???.tgz and extract it in the solvers/osi/ folder.
The version that the Makefile presumes is Osi-0.105.2. If you download a different version the line `OSIVER = 0.105.2` will have to be changed.

#### OsiClp
Download Clp from: [http://www.coin-or.org/download/source/Clp/][osiclp]

The version that the Makefile presumes is Clp-1.14.5. If you download a different version the line `CLPVER = 1.14.5` will have to be changed.

#### OsiCbc
Download Cbc from: [http://www.coin-or.org/download/source/Cbc/][osicbc]

The version that the Makefile presumes is Cbc-2.7.5. If you download a different version the line `CBCVER = 2.7.5` will have to be changed.

CBC Supports several different lp solver backends(_e.g. glpk_) but currently dynamically links Clp.
There will be some work done on this, but as it stands you need to have Clp installed to the system for CBC to work.

#### OsiVol
Download Vol from: [http://www.coin-or.org/download/source/Vol/][osivol]

The version that the Makefile presumes is Vol-1.3.2. If you download a different version the line `VOLVER = 1.3.2` will have to be changed.

#### OsiDylp
Download Dylp from: [http://www.coin-or.org/download/source/Dylp/][osidylp]

The version that the Makefile presumes is DyLP-1.8.2. If you download a different version the line `DYLPVER = 1.8.2` will have to be changed.

#### OsiSpx (soplex)
Download soplex from: [http://soplex.zib.de/download.shtml][soplex]

The version that the Makefile presumes is soplex-1.6.0. If you download a different version the line `SPXVER = 1.6.0` will have to be changed.

#### OsiGlpk
Download Glpk from: [http://ftp.gnu.org/gnu/glpk/][glpk]

The version that the Makefile presumes is glpk-4.47. If you download a different version the line `GLPKVER = 4.47` will have to be changed.

[njhome]: http://numberjack.ucc.ie
[osi]: http://www.coin-or.org/download/source/Osi/
[osiclp]: http://www.coin-or.org/download/source/Clp/
[osicbc]: http://www.coin-or.org/download/source/Cbc/
[osivol]: http://www.coin-or.org/download/source/Vol/
[osidylp]: http://www.coin-or.org/download/source/DyLP/
[soplex]: http://soplex.zib.de/download.shtml
[glpk]: http://ftp.gnu.org/gnu/glpk/
