Contributing to Numberjack
==========================

This document outlines various aspects for contributing to Numberjack, such as requirements to
compile locally, coding standards, testing, and building the documentation.


.. toctree::


Requirements
------------

Numberjack itself is mostly written in Python, with some of the interfaces to solvers writting in
C++. In order to build Numberjack and some of the underlying solver that are included, the following
requirements will need to::

    sudo apt-get install build-essential python3-dev swig libxml2-dev zlib1g-dev libgmp-dev



Code Standards
--------------

Although much of the older code does not follow these guidelines, please try to follow these for
new additions. This makes it more consistent going forward, and easier for others to follow.

Python code should follow the `PEP8`_ guidlines, with the exception of line-lengths being limited to
100 characters instead of 80, where possible.

.. _PEP8: https://www.python.org/dev/peps/pep-0008/


Testing
-------


Unit-tests
^^^^^^^^^^

Nose is used to run the unit-tests. To install: `pip install nose`. To run the unit-tests, from
the Numberjack root folder, run `nosetests tests`, where `tests` is the folder containint the
unit-tests.

Multiple environments
^^^^^^^^^^^^^^^^^^^^^

Tox is used to test the build, install, and unit-tests on multiple python environments. To install
it: `pip install tox`. The versions of python specified in `tox.ini` should be installed on the
system. Then, to run the tests: `tox`. This will launch a new virtual environment for each python
version, build a source distribution, install it, and run the unit-tests. The output of the test
command will be printed to screen as well as any build errors.


Documentation
-------------

To build the documentation locally, you will need to install Sphinx: `pip install sphinx`
Then `cd` into the `doc` directory and run `make` with the format you want to produce,
such as `make html`. This creates the HTML version of the documentation in `doc/build/html`.

This should be copied to the numberjack.ucc.ie server and placed in the relevant directory.



Building a new distribution
---------------------------

Make sure your working directory is clean and up to date with the github repository. This avoids
leaking unwanted files into the distribution, because all files in certain folders, like `examples`,
will be included in the distribution.

Building a new source distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To produce a zipped file in the `dist` folder with the source distribution::

    python setup.py sdist

To upload the source distribution to the Python Package Index (PyPI)::

    twine upload dist/Numberjack-?.?.?.tar.gz  # replacing the ? with the version number
