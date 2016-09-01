# toulbar2
## An exact solver for cost function networks

## What is toulbar2?

toulbar2 is an  open-source C++ solver for cost  function networks. It
solves  various combinatorial  optimization problems.  The constraints
and objective function  are factorized in local  functions on discrete
variables. Each  function returns a  cost (a finite  positive integer)
for any  assignment of its  variables. Constraints are  represented as
functions with costs in {0,k} where k is a large integer representing
forbidden assignments.  toulbar2 looks for a  non-forbidden assignment
of all variables  that minimizes the sum of all  functions. 

toulbar2 won  several competitions on deterministic  and probabilistic
graphical models:

* Max-CSP 2008 Competition [CPAI08][cpai08] (winner on 2-ARY-EXT and N-ARY-EXT)
* Probabilistic Inference Evaluation [UAI 2008][uai2008] (winner on several MPE tasks, inra entries)
* 2010 UAI APPROXIMATE INFERENCE CHALLENGE [UAI 2010][uai2010] (winner on 1200-second MPE task)
* The Probabilistic Inference Challenge [PIC 2011][pic2011] (second place by ficolofo on 1-hour MAP task)
* UAI 2014 Inference Competition [UAI 2014][uai2014] (winner on all MAP task categories, see Proteus, Robin, and IncTb entries)

[cpai08]: http://www.cril.univ-artois.fr/CPAI08/
[uai2008]: http://graphmod.ics.uci.edu/uai08/Evaluation/Report
[uai2010]: http://www.cs.huji.ac.il/project/UAI10/summary.php
[pic2011]: http://www.cs.huji.ac.il/project/PASCAL/board.php
[uai2014]: http://www.hlt.utdallas.edu/~vgogate/uai14-competition/leaders.html 


## Download

http://mulcyber.toulouse.inra.fr/projects/toulbar2/

Latest src/debian/windows x86_64 releases:
* 2016: [src0_9_8]/[deb0_9_8]/[win0_9_8] with hybrid best-first search and more soft global cost functions
* 2015: [src0_9_7]/[deb0_9_7]/[win0_9_7] with local search [INCOP][incop] solver after preprocessing
* 2014: [src0_9_6]/[deb0_9_6]/[win0_9_6] with extra variable ordering heuristics and dominance pruning rules
* 2012: [src0_9_5]/[deb0_9_5]/[win0_9_5] with soft global decomposable cost functions
* 2011: [src0_9_4]/[deb0_9_4]/[win0_9_4] with more preprocessing techniques
* 2010: [src0_9_3]/[deb0_9_3]/[win0_9_3] with soft global cost functions

[src0_9_8]: https://mulcyber.toulouse.inra.fr/frs/download.php/1455/toulbar2.0.9.8.0-Release-sources.tar.gz
[src0_9_7]: https://mulcyber.toulouse.inra.fr/frs/download.php/1380/toulbar2.0.9.7.0-Release-sources.tar.gz
[src0_9_6]: https://mulcyber.toulouse.inra.fr/frs/download.php/1292/toulbar2.0.9.6.0-Release-sources.tar.gz
[src0_9_5]: https://mulcyber.toulouse.inra.fr/frs/download.php/1142/toulbar2.0.9.5.0-Release-sources.tar.gz
[src0_9_4]: https://mulcyber.toulouse.inra.fr/frs/download.php/1019/toulbar2.0.9.4.0-Release-sources.tar.gz
[src0_9_3]: https://mulcyber.toulouse.inra.fr/frs/download.php/975/toulbar2.0.9.3.0-Release-sources.tar.gz

[deb0_9_3]: https://mulcyber.toulouse.inra.fr/frs/download.php/964/toulbar2.0.9.3.0-Release-i686.deb
[deb0_9_4]: https://mulcyber.toulouse.inra.fr/frs/download.php/1008/toulbar2.0.9.4.0-Release-i686.deb
[deb0_9_5]: https://mulcyber.toulouse.inra.fr/frs/download.php/1134/toulbar2.0.9.5.0-Release-x86_64.deb
[deb0_9_6]: https://mulcyber.toulouse.inra.fr/frs/download.php/1281/toulbar2.0.9.6.0-Release-i686.deb
[deb0_9_7]: https://mulcyber.toulouse.inra.fr/frs/download.php/1371/toulbar2.0.9.7.0-Release-x86_64.deb
[deb0_9_8]: https://mulcyber.toulouse.inra.fr/frs/download.php/1448/toulbar2.0.9.8.0-Release-x86_64.deb

[win0_9_3]: https://mulcyber.toulouse.inra.fr/frs/download.php/962/toulbar2.0.9.3.0-Release-i686.exe
[win0_9_4]: https://mulcyber.toulouse.inra.fr/frs/download.php/1006/toulbar2.0.9.4.0-Release-i686.exe
[win0_9_5]: https://mulcyber.toulouse.inra.fr/frs/download.php/1129/toulbar2.0.9.5.0-Release-i686.exe
[win0_9_6]: https://mulcyber.toulouse.inra.fr/frs/download.php/1279/toulbar2.0.9.6.0-Release-i686.exe
[win0_9_7]: https://mulcyber.toulouse.inra.fr/frs/download.php/1374/toulbar2.0.9.7.0-Release-x86_64.exe
[win0_9_8]: https://mulcyber.toulouse.inra.fr/frs/download.php/1446/toulbar2.0.9.8.0-Release-x86_64.exe


## Installation

Library needed:
* libgmp-dev
* libboost-dev
* libboost-graph-dev

Optional libraries:
* libxml2-dev

GNU C++ Symbols to be defined if using Linux Eclipse/CDT IDE (no value needed):
* LINUX
* LONGLONG_COST
* WIDE_STRING
* LONGDOUBLE_PROB
* NARYCHAR
* WCSPFORMATONLY

Commands for compiling toulbar2 on Linux in directory toulbar2/src without cmake:

    bash
    cd src
    echo '#define Toulbar_VERSION "0.9.8"' > ToulbarVersion.hpp
    g++ -o toulbar2 -I. tb2*.cpp incop/*.cpp ToulbarVersion.cpp -O3 -std=c++11 -DNDEBUG -DLINUX \
     -DLONGLONG_COST -DWIDE_STRING -DLONGDOUBLE_PROB -DNARYCHAR -DWCSPFORMATONLY -lgmp -static

## Authors

toulbar2 was originally developped by Toulouse (INRA MIAT) and Barcelona (UPC, IIIA-CSIC) teams, hence the name of the solver. 

Additional contributions by:
* The Chinese University of Hong Kong and Caen University, France (GREYC) for global cost functions
* Marseille University, France (LSIS) for tree decomposition heuristics
* Ecole des Ponts ParisTech, France (CERMICS/LIGM) for [INCOP][incop] local search solver
* University College Cork, Ireland (Insight) for a Python interface in [NumberJack][numberjack] and a portfolio dedicated to UAI graphical models [Proteus][proteus]
* Artois University, France (CRIL) for an XCSP 2.1 format reader of CSP and WCSP instances

[incop]: http://imagine.enpc.fr/~neveub/incop/incoppresentation.html
[numberjack]: http://numberjack.ucc.ie/
[proteus]: https://github.com/9thbit/uai-proteus


## Citing

Please use one of the following references for citing toulbar2:

 Multi-Language Evaluation of Exact Solvers in Graphical Model Discrete Optimization
 Barry Hurley, Barry O'Sullivan, David Allouche, George Katsirelos, Thomas Schiex, Matthias Zytnicki, Simon de Givry
 Constraints, 21(3):413-434, 2016

 Tractability-preserving Transformations of Global Cost Functions
 David Allouche, Christian Bessiere, Patrice Boizumault, Simon de Givry, Patricia Gutierrez, Jimmy HM. Lee, Ka Lun Leung, Samir Loudni, Jean-Philippe Métivier, Thomas Schiex, Yi Wu
 Artificial Intelligence, 238:166-189, 2016

 Soft arc consistency revisited
 Martin Cooper, Simon de Givry, Marti Sanchez, Thomas Schiex, Matthias Zytnicki, and Thomas Werner
 Artificial Intelligence, 174(7-8):449-478, 2010 


##  What are the algorithms inside toulbar2?

* Soft arc consistencies (NC, AC, DAC, FDAC)
 In the quest of the best form of local consistency for Weighted CSP
 J. Larrosa & T. Schiex
 In Proc. of IJCAI-03. Acapulco, Mexico, 2003

* Soft existential arc consistency (EDAC)
 Existential arc consistency: Getting closer to full arc consistency in weighted csps
 S. de Givry, M. Zytnicki, F. Heras, and J. Larrosa
 In Proc. of IJCAI-05, Edinburgh, Scotland, 2005

* Depth-first Branch and Bound exploiting a tree decomposition (BTD)
 Exploiting Tree Decomposition and Soft Local Consistency in Weighted CSP
 S. de Givry, T. Schiex, and G. Verfaillie
 In Proc. of AAAI-06, Boston, MA, 2006 

* Virtual arc consistency (VAC)
 Virtual arc consistency for weighted csp
 M. Cooper, S. de Givry, M. Sanchez, T. Schiex, and M. Zytnicki
 In Proc. of AAAI-08, Chicago, IL, 2008

* Soft generalized arc consistencies (GAC, FDGAC)
 Towards Efficient Consistency Enforcement for Global Constraints in Weighted Constraint Satisfaction
 J. H. M. Lee and K. L. Leung
 In Proc. of IJCAI-09, Los Angeles, USA, 2010

* Russian doll search exploiting a tree decomposition (RDS-BTD)
 Russian doll search with tree decomposition
 M Sanchez, D Allouche, S de Givry, and T Schiex
 In Proc. of IJCAI'09, Pasadena (CA), USA, 2009

* Soft bounds arc consistency (BAC)
 Bounds Arc Consistency for Weighted CSPs
 M. Zytnicki, C. Gaspin, S. de Givry, and T. Schiex
 Journal of Artificial Intelligence Research, 35:593-621, 2009

* Counting solutions in satisfaction (#BTD, Approx_#BTD)
 Exploiting problem structure for solution counting
 A. Favier, S. de Givry, and P. Jégou
 In Proc. of CP-09, Lisbon, Portugal, 2009

* Soft existential generalized arc consistency (EDGAC)
 A Stronger Consistency for Soft Global Constraints in Weighted Constraint Satisfaction
 J. H. M. Lee and K. L. Leung
 In Proc. of AAAI-10, Boston, MA, 2010 

* Preprocessing techniques (combines variable elimination and cost function decomposition)
 Pairwise decomposition for combinatorial optimization in graphical models
 A Favier, S de Givry, A Legarra, and T Schiex
 In Proc. of IJCAI-11, Barcelona, Spain, 2011

* Decomposable global cost functions (wregular, wamong, wsum) 
 Decomposing global cost functions
 D Allouche, C Bessiere, P Boizumault, S de Givry, P Gutierrez, S Loudni, JP Métivier, and T Schiex
 In Proc. of AAAI-12, Toronto, Canada, 2012

* Pruning by dominance (DEE)
 Dead-End Elimination for Weighted CSP
 S de Givry, S Prestwich, and B O'Sullivan
 In Proc. of CP-13, pages 263-272, Uppsala, Sweden, 2013

* Hybrid best-first search exploiting a tree decomposition (HBFS)
 Anytime Hybrid Best-First Search with Tree Decomposition for Weighted CSP
 D Allouche, S de Givry, G Katsirelos, T Schiex, and M Zytnicki
 In Proc. of CP-15, Cork, Ireland, 2015 


Copyright (C) 2006-2016, INRA.
toulbar2 is currently maintained by Simon de Givry, INRA - MIAT, Toulouse, France (simon.degivry@toulouse.inra.fr)

