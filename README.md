# Numberjack
## A python constraint programming platform 

## What is Numberjack?

From the Numberjack [homepage][njhome]

Numberjack is a modelling package written in Python for constraint programming. Python benefits from a large and active programming community, Numberjack is therefore a perfect tool to embed CP technology into larger applications. It is designed to support a number of underlying C/C++ solvers as egg files, that is, seamlessly and efficiently. Currently, there are four available back-ends: a MIP solver (SCIP), two SAT solvers (MiniSat) and (Walksat) and CP solver (Mistral).

* Numberjack offers a high level constraint programming language
* Numberjack directly benefits from python's features and modules
* Numberjack uses efficient underlying C/C++ solvers, without the compilation hassle!

## Building the solvers

### Selecting solvers

At the moment there is no mechanism to select which solvers get compiled so deleting the solvers you dont want or moving them out of the solvers folder is needed.

### Prerequisites

* python 2 (`python` should link to python2, arch linux uses python3 by default which makes things interesting, a python 2 virtualenv helps, note python-config must point to python2-config)
* swig
* libxml2-dev
* zlib1g-dev
* python-dev

### Mistral, Minisat, Walksat

The sourcecode for these solvers is included in the source tree.
Running make should suffice

### SCIP
Due to licencing issues Numberjack cannot have the sources for SCIP included in the source tree.
The SCIP source code is available under an Academic Licence for free from:
[http://scip.zib.de/download.shtml][scipopt]
Please download version 3.0.0.

1. Download the source code of SCIP Optimization Suite v3.0.0, unarchive it, however, do not install nor compile it.

```bash
    tar zxf scipoptsuite-3.0.0.tgz
```

2. In the scipoptsuite directory, unarchive the SCIP and SoPlex archives.

```bash
    cd scipoptsuite-3.0.0/
    tar zxf scip-3.0.0.tgz
    tar zxf soplex-1.7.0.tgz
```

3. Set an environment variable ZIBPATH to the scipoptsuite directory, for example:

```bash
    export ZIBPATH=path_to/scipoptsuite-3.0.0
```

### Osi Solvers
Some of the solvers cannot be included in the source tree. I've left out all of them to be safe.
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
[scipopt]: http://scip.zib.de/download.shtml
[osi]: http://www.coin-or.org/download/source/Osi/
[osiclp]: http://www.coin-or.org/download/source/Clp/
[osicbc]: http://www.coin-or.org/download/source/Cbc/
[osivol]: http://www.coin-or.org/download/source/Vol/
[osidylp]: http://www.coin-or.org/download/source/DyLP/
[soplex]: http://soplex.zib.de/download.shtml
[glpk]: http://ftp.gnu.org/gnu/glpk/
