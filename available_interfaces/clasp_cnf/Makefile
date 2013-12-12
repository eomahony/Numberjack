include ../../make.project

MAINDIR = $(shell pwd)/../..
CURDIR ?= .
PYTHON = $(CURDIR)/python

PYFILE      =   $(PYTHON)/ClaspCNF.py

wrapper:

clean:
	@rm -rf *~ $(PYTHON)/*.pyc

clean_swig:

install_python: wrapper
	cd $(PYTHON); python $(MAINDIR)/tools/setup.py install

externalsolver_local: wrapper
