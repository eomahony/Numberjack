#! /usr/bin/env python

'''
  Numberjack is a constraint satisfaction and optimisation library
  Copyright (C) 2009-2013 Cork Constraint Computation Center, UCC
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.
  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The authors can be contacted electronically at 
  numberjack.support@gmail.com
'''


from distutils.core import setup
setup(name='Numberjack',
      version='1.0',
      author="Numberjack Team",
      author_email="numberjack.support@gmail.com",
      url="http://numberjack.ucc.ie/",
      py_modules=['Numberjack', 'MIPParser', 'ExternalSolver', 'XCSP', 'XCSPOut', 'Decomp'],
      license="License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
      description="A Python platform for combinatorial optimization.",
      long_description="""Numberjack is a modelling package written in Python
for constraint programming and combinatorial optimization. Python benefits from
a large and active programming community, Numberjack is therefore a perfect tool
to embed CP technology into larger applications. It is designed to support a
number of efficient underlying C/C++ solvers seamlessly and efficiently.""",
      )
