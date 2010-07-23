'''
  Numberjack is a constraint satisfaction and optimisation library
  Copyright (C) 2009 Emmanuel Hebrard, Eoin O'Mahony
  
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


#! /usr/bin/env python

from distutils.core import setup
setup(name='Numberjack',
      version='1.0',
      py_modules=['Numberjack'],
      )
setup(name='Decomp',
      version='1.0',
      py_modules=['Decomp'],
      )

#setup(name="ThreadSolver",
#      version='1.0',
#      py_modules=['ThreadSolver'])
#setup(name='NumberjackSolver',
#      version='1.0',
#      package_dir={'NumberjackSolver': 'NumberjackSolver'},
#      packages=['NumberjackSolver']
#      )
