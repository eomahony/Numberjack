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

from __future__ import print_function
from distutils.core import setup
from distutils.extension import Extension
from distutils.command.build_ext import build_ext as _build_ext
from distutils.errors import CCompilerError
from distutils.spawn import find_executable
import subprocess as sp
import sys
import os


EXTRA_COMPILE_ARGS = [
    '-O3',
    # ignore warnings about swig code
    '-Wno-shadow', '-Wno-unused-label',
]
EXTRA_LINK_ARGS = []
extensions = []  # List of solver interfaces (C extensions) to be built.
disabled_extensions = []  # The names of solver interfaces which have been
    # disabled because the solver was not found on the system.


if sys.platform == 'darwin':
    EXTRA_COMPILE_ARGS.extend([
        '-stdlib=libstdc++', '-Wno-shorten-64-to-32',
        '-arch', 'x86_64',  # force 64bit only builds on Mac

        # ignore warning about swig code
        '-Wno-self-assign', '-Wno-unused-const-variable',
    ])
    EXTRA_LINK_ARGS.extend(['-arch', 'x86_64', '-stdlib=libstdc++'])
elif sys.platform.startswith('linux'):
    EXTRA_COMPILE_ARGS.extend([
        '-Wno-narrowing',

        # ignore warning about swig code
        '-Wno-unused-but-set-variable',
    ])


class njbuild_ext(_build_ext):
    # Allow for C extension compilation to fail.

    def __init__(self, *args, **kwargs):
        self.builtsolvernames = []
        self.failedsolvernames = []
        _build_ext.__init__(self, *args, **kwargs)

    def run(self):
        _build_ext.run(self)

        if disabled_extensions:
            print("The following solvers could not be located "
                  "on your system so their interface has been disabled:",
                  ", ".join(disabled_extensions), file=sys.stderr)

        if self.failedsolvernames:
            print("Failed to build the following interfaces "
                  "(details are above):", ", ".join(self.failedsolvernames), file=sys.stderr)

        if self.builtsolvernames:
            print("Successfully built solver interfaces for",
                  ", ".join(self.builtsolvernames))

    def build_extension(self, ext):
        try:
            _build_ext.build_extension(self, ext)

            # Record the names of extensions which were successful, trim '_'
            self.builtsolvernames.append(ext.name[1:])

        except CCompilerError:
            self.failedsolvernames.append(ext.name[1:])


# ------------------------------ Helper Methods ------------------------------


def xml2config(option, path="xml2-config"):
    import shlex
    if not find_executable(path):
        raise RuntimeError("Could not find %s" % path)
    cmd = '%s %s' % (path, option)
    p = sp.Popen(cmd, stdout=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    return shlex.split(stdout.decode("utf-8").strip())


CPLEX, GUROBI, SCIP = "CPLEX", "Gurobi", "SCIP"


def get_solver_home(solvername):
    # Helps to find the path to a particular third party solver

    if solvername == CPLEX:
        # Try for environmental variable first
        env_path = os.getenv('CPLEXDIR')
        if env_path and len(env_path.strip()) > 0:
            return env_path

        # Try to find the cplex binary in the PATH
        ex_path = find_executable('cplex')
        if ex_path:
            ex_path = os.path.realpath(ex_path)  # Expand symbolic links if any
            ex_dir = os.path.dirname(ex_path)  # Path to the bin directory
            return os.path.abspath(os.path.join(ex_dir, os.pardir, os.pardir))

    elif solvername == GUROBI:
        # Try for environmental variable first
        env_path = os.getenv('GUROBI_HOME')
        if env_path and len(env_path.strip()) > 0:
            return env_path

        # Try to find the gurobi_cl binary in the PATH
        ex_path = find_executable('gurobi_cl')
        if ex_path:
            ex_path = os.path.realpath(ex_path)  # Expand symbolic links if any
            ex_dir = os.path.dirname(ex_path)  # Path to the bin directory
            return os.path.abspath(os.path.join(ex_dir, os.pardir))

    elif solvername == SCIP:
        # Try for environmental variable first
        env_path = os.getenv('ZIBPATH')
        if env_path and len(env_path.strip()) > 0:
            return env_path

        # Otherwise, if the tgz has been placed in this directory, extract it
        # and return that folder's path
        scipoptfolder = None
        for basename in os.listdir(os.getcwd()):
            if os.path.isfile(basename) and \
                    basename.startswith("scipoptsuite-") and \
                    basename.endswith(".tgz"):
                scipoptfolder = os.path.join(os.getcwd(), basename[:-4])
                if not os.path.isdir(scipoptfolder):
                    sp.call('tar zxf %s' % basename,
                            cwd=os.getcwd(), shell=True)
                break
        return scipoptfolder

    else:
        raise RuntimeError("Error unknown solver name '%s'" % solvername)

    return None


# ------------------------------ Extensions ------------------------------


def listextfiles(folder, extension=".cpp"):
    ret = []
    for bn in os.listdir(folder):
        if os.path.splitext(bn)[1] == extension:
            ret.append(os.path.join(folder, bn))
    return ret


mistralsrc = 'Numberjack/solvers/Mistral/mistral/lib/src'
mistral = Extension(
    '_Mistral',
    sources=[
        'Numberjack/solvers/Mistral.i',
        'Numberjack/solvers/Mistral/Mistral.cpp',
    ] + listextfiles(mistralsrc),
    swig_opts=[
        '-modern', '-c++',
        '-INumberjack/solvers/Mistral',
    ],
    include_dirs=[
        'Numberjack/solvers/Mistral',
        'Numberjack/solvers/Mistral/mistral/include',
    ],
    libraries=['m'],
    language='c++',
    define_macros=[('_UNIX', None)],
    extra_compile_args=EXTRA_COMPILE_ARGS + xml2config('--cflags') +
    ['-fPIC', '-Wno-unused-label', '-fexceptions', '-Wno-sequence-point',
     '-Wno-overloaded-virtual', '-Wno-narrowing', '-Wno-unused-variable',
     '-Wno-unused-but-set-variable'],
    extra_link_args=EXTRA_LINK_ARGS + xml2config('--libs'),
)
extensions.append(mistral)


mistral2src = 'Numberjack/solvers/Mistral2/mistral/src/lib'
mistral2 = Extension(
    '_Mistral2',
    sources=[
        'Numberjack/solvers/Mistral2.i',
        'Numberjack/solvers/Mistral2/Mistral2.cpp',
    ] + listextfiles(mistral2src),
    swig_opts=[
        '-modern', '-c++',
        '-INumberjack/solvers/Mistral2',
    ],
    include_dirs=[
        'Numberjack/solvers/Mistral2',
        'Numberjack/solvers/Mistral2/mistral/src/include',
        'Numberjack/solvers/Mistral2/mistral/tools/tclap/include'
    ],
    libraries=['m'],
    language='c++',
    # define_macros=[('_UNIX', None)],
    extra_compile_args=EXTRA_COMPILE_ARGS +
    ['-fPIC', '-Wno-unused-label', '-fexceptions', '-Wno-overloaded-virtual',
     '-Wno-unused-variable', '-Wno-parentheses', '-Wno-reorder',
     '-Wno-delete-non-virtual-dtor'],
    extra_link_args=EXTRA_LINK_ARGS,
)
extensions.append(mistral2)


toulbar2src = 'Numberjack/solvers/Toulbar2/lib/src'
toulbar2incopsrc = 'Numberjack/solvers/Toulbar2/lib/src/incop'
toulbar2 = Extension(
    '_Toulbar2',
    sources=[
        'Numberjack/solvers/Toulbar2.i',
        'Numberjack/solvers/Toulbar2/Toulbar2.cpp',
    ] + listextfiles(toulbar2src) + listextfiles(toulbar2incopsrc),
    swig_opts=[
        '-modern', '-c++',
        '-INumberjack/solvers/Toulbar2',
    ],
    include_dirs=[
        'Numberjack/solvers/Toulbar2',
        'Numberjack/solvers/Toulbar2/include',
        'Numberjack/solvers/Toulbar2/include/incop',
    ],
    libraries=['gmp'],
    language='c++',
    define_macros=[
        ('NUMBERJACK', None),
        # ('BOOST', None), # requires libboost-graph-dev installed
        ('NDEBUG', None),
        ('LINUX', None),
        ('LONGLONG_COST', None),
        ('WIDE_STRING', None),
        ('LONGDOUBLE_PROB', None),
        ('NARYCHAR', None),
    ],
    extra_compile_args=EXTRA_COMPILE_ARGS + [
        '-Wno-overloaded-virtual', '-Wno-mismatched-tags', '-Wno-unused-private-field',
    ],
    extra_link_args=EXTRA_LINK_ARGS,
)
extensions.append(toulbar2)


mip = Extension(
    '_MipWrapper',
    sources=[
        'Numberjack/solvers/MipWrapper.i',
        'Numberjack/solvers/MipWrapper/MipWrapper.cpp',
    ],
    swig_opts=[
        '-modern', '-c++',
        '-INumberjack/solvers/MipWrapper',
    ],
    include_dirs=['Numberjack/solvers/MipWrapper'],
    language='c++',
    extra_compile_args=EXTRA_COMPILE_ARGS,
    extra_link_args=EXTRA_LINK_ARGS,
)
extensions.append(mip)


sat = Extension(
    '_SatWrapper',
    sources=[
        'Numberjack/solvers/SatWrapper.i',
        'Numberjack/solvers/SatWrapper/SatWrapper.cpp',
    ],
    swig_opts=[
        '-modern', '-c++',
        '-INumberjack/solvers/SatWrapper',
        '-INumberjack/solvers/MiniSat/minisat_src/core',
        '-INumberjack/solvers/MiniSat/minisat_src/mtl',
    ],
    include_dirs=[
        'Numberjack/solvers/SatWrapper/',
        'Numberjack/solvers/MiniSat/minisat_src/core',
        'Numberjack/solvers/MiniSat/minisat_src/mtl'
    ],
    language='c++',
    extra_compile_args=EXTRA_COMPILE_ARGS,
    extra_link_args=EXTRA_LINK_ARGS,
)
extensions.append(sat)


minisat = Extension(
    '_MiniSat',
    sources=[
        'Numberjack/solvers/MiniSat.i',
        'Numberjack/solvers/MiniSat/MiniSat.cpp',
        'Numberjack/solvers/SatWrapper/SatWrapper.cpp',
        'Numberjack/solvers/MiniSat/SimpSolver.cpp',
        'Numberjack/solvers/MiniSat/minisat_src/core/Solver.C',
    ],
    swig_opts=[
        '-modern', '-c++',
        '-INumberjack/solvers/MiniSat',
        '-INumberjack/solvers/SatWrapper',
        '-INumberjack/solvers/MiniSat/minisat_src/core',
        '-INumberjack/solvers/MiniSat/minisat_src/mtl',
    ],
    include_dirs=[
        'Numberjack/solvers/MiniSat',
        'Numberjack/solvers/SatWrapper',
        'Numberjack/solvers/MiniSat/minisat_src/core',
        'Numberjack/solvers/MiniSat/minisat_src/mtl'
    ],
    language='c++',
    extra_compile_args=EXTRA_COMPILE_ARGS,
    extra_link_args=EXTRA_LINK_ARGS,
)
extensions.append(minisat)


walksat = Extension(
    '_Walksat',
    sources=[
        'Numberjack/solvers/Walksat.i',
        'Numberjack/solvers/Walksat/Walksat.cpp',
        'Numberjack/solvers/Walksat/walksat_src/cpp_walksat.cpp',
        'Numberjack/solvers/SatWrapper/SatWrapper.cpp',
    ],
    swig_opts=[
        '-modern', '-c++',
        '-INumberjack/solvers/Walksat',
        '-INumberjack/solvers/SatWrapper',
        '-INumberjack/solvers/Walksat/walksat_src',
        '-INumberjack/solvers/MiniSat/minisat_src/core',
        '-INumberjack/solvers/MiniSat/minisat_src/mtl',
    ],
    include_dirs=[
        'Numberjack/solvers/Walksat',
        'Numberjack/solvers/SatWrapper',
        'Numberjack/solvers/Walksat/walksat_src',
        'Numberjack/solvers/MiniSat/minisat_src/core',
        'Numberjack/solvers/MiniSat/minisat_src/mtl'
    ],
    language='c++',
    extra_compile_args=EXTRA_COMPILE_ARGS +
    ['-ffloat-store', '-Wno-format', '-Wno-unused-variable'],
    extra_link_args=EXTRA_LINK_ARGS,
)
extensions.append(walksat)


cplexhome = get_solver_home(CPLEX)
if cplexhome:
    concertdir = os.path.abspath(os.path.join(cplexhome, os.pardir, "concert"))

    def get_cplex_includes():
        cplexincdir = os.path.join(cplexhome, "include")
        concertincdir = os.path.join(concertdir, "include")
        return [cplexincdir, concertincdir]

    def get_cplex_lib_dirs():
        """
        CPLEX lib structure: cplexdir/lib/ARCHITECTURE/BUILDTYPE/
        This function will walk the cplexdir/lib folder to find the
        folders containing static libraries.
        """

        def getlibdirs(libfolder):
            for dirpath, dirnames, filenames in os.walk(libfolder):
                for f in filenames:
                    if f.endswith(".a"):
                        return str(dirpath)
            raise RuntimeError(
                "Error could not find the lib folder in %s" % cplexhome)

        cplexlibfolder = os.path.join(cplexhome, "lib")
        concertlibfolder = os.path.join(concertdir, "lib")
        return [getlibdirs(cplexlibfolder), getlibdirs(concertlibfolder)]

    cplex = Extension(
        '_CPLEX',
        sources=[
            'Numberjack/solvers/CPLEX.i',
            'Numberjack/solvers/CPLEX/CPLEX.cpp',
            'Numberjack/solvers/MipWrapper/MipWrapper.cpp',
        ],
        swig_opts=[
            '-modern', '-c++',
            '-INumberjack/solvers/CPLEX',
            '-INumberjack/solvers/MipWrapper',
        ],
        include_dirs=[
            'Numberjack/solvers/CPLEX',
            'Numberjack/solvers/MipWrapper',
        ] + get_cplex_includes(),
        library_dirs=get_cplex_lib_dirs(),
        language='c++',
        define_macros=[('_UNIX', None), ('NDEBUG', None), ('IL_STD', None)],
        extra_compile_args=EXTRA_COMPILE_ARGS +
        ['-O', '-fPIC', '-fexceptions'] +
        ["-fno-strict-aliasing"] if sys.platform.startswith('linux') else [],
        libraries=['concert', 'ilocplex', 'cplex', 'm', 'pthread'],
        extra_link_args=EXTRA_LINK_ARGS,
    )
    extensions.append(cplex)
else:
    disabled_extensions.append(CPLEX)


gurobihome = get_solver_home(GUROBI)
if gurobihome:
    gurobiincdir = os.path.join(gurobihome, "include")
    gurobilibdir = os.path.join(gurobihome, "lib")

    def get_gurobi_libs():

        def get_gurobi_libname():
            import re
            libre = re.compile("lib(?P<libname>gurobi\d+)\.so")
            for dirpath, dirnames, filenames in os.walk(gurobilibdir):
                for f in filenames:
                    match = libre.match(f)
                    if match:
                        return match.groupdict()["libname"]
            raise RuntimeError(
                "Error could not find the Gurobi library in '%s'" %
                gurobilibdir)
        return ['gurobi_c++', get_gurobi_libname()]

    gurobi = Extension(
        '_Gurobi',
        sources=[
            'Numberjack/solvers/Gurobi.i',
            'Numberjack/solvers/Gurobi/Gurobi.cpp',
            'Numberjack/solvers/MipWrapper/MipWrapper.cpp',
        ],
        swig_opts=[
            '-modern', '-c++',
            '-INumberjack/solvers/Gurobi',
            '-INumberjack/solvers/MipWrapper',
        ],
        include_dirs=[
            'Numberjack/solvers/Gurobi',
            'Numberjack/solvers/MipWrapper',
            gurobiincdir,
        ],
        library_dirs=[gurobilibdir],
        libraries=get_gurobi_libs(),
        extra_compile_args=EXTRA_COMPILE_ARGS +
        ['-fPIC', '-fexceptions'],
        extra_link_args=EXTRA_LINK_ARGS,
        language='c++',
    )
    extensions.append(gurobi)
else:
    disabled_extensions.append(GUROBI)


scipopthome = get_solver_home(SCIP)
if scipopthome:
    scipoptlibfolder = os.path.join(scipopthome, "lib")

    def get_scip_libs():
        import re
        libre = re.compile("lib(?P<libname>scipopt.*)\.a")
        for dirpath, dirnames, filenames in os.walk(scipoptlibfolder):
            for f in filenames:
                match = libre.match(f)
                if match:
                    return [match.groupdict()["libname"]]
        raise RuntimeError(
            "Error could not find the SCIP library in '%s'" %
            scipoptlibfolder)

    def compile_scip_ifneeded():
        # Check if the static lib already exists and if not compile it.
        # There may be a better place to do this other, e.g. only for build_ext
        try:
            get_scip_libs()
            return
        except RuntimeError:
            makecmd = "make scipoptlib ZIMPL=false ZLIB=false READLINE=false" \
                " GAMS=false GMP=false LEGACY=true SPX_LEGACY=true"

            print("Compiling SCIP library...")
            returncode = sp.call(makecmd, cwd=scipopthome, shell=True)
            if returncode != 0:
                sys.exit(1)

    def get_scipcomponent_folder(comp):
        for basename in os.listdir(scipopthome):
            fullpath = os.path.join(scipopthome, basename)
            if os.path.isdir(fullpath) and basename.startswith(comp):
                return fullpath

        raise RuntimeError("Error could not find the SCIP %s folder." % comp)

    compile_scip_ifneeded()

    scipfolder = get_scipcomponent_folder("scip")
    scipincfolder = os.path.join(scipfolder, "src")

    scip = Extension(
        '_SCIP',
        sources=[
            'Numberjack/solvers/SCIP.i',
            'Numberjack/solvers/SCIP/SCIP.cpp',
            'Numberjack/solvers/MipWrapper/MipWrapper.cpp',
        ],
        swig_opts=[
            '-modern', '-c++',
            '-INumberjack/solvers/SCIP',
            '-INumberjack/solvers/MipWrapper',
        ],
        include_dirs=[
            'Numberjack/solvers/SCIP',
            'Numberjack/solvers/MipWrapper',
            scipincfolder,
        ],
        library_dirs=[scipoptlibfolder],
        libraries=get_scip_libs(),
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        language='c++',
    )
    extensions.append(scip)
else:
    disabled_extensions.append(SCIP)


# ------------------------------ End Extensions ------------------------------

# Possibly compile only a subset of the solvers. This can be specified as a CSV
# list on the command line or by passing -solver multiple times, e.g:
# python setup.py build -solver Mistral,Mistral2 -solver SCIP
allsolvers = dict((e.name[1:], e) for e in extensions)
SOLVERARG = "-solver"
solversubsetnames = set()
while SOLVERARG in sys.argv:
    pos = sys.argv.index(SOLVERARG)
    if pos + 1 >= len(sys.argv):
        print("Error: you should specify a solver or list of "
              "solvers to compile with the %s option." % SOLVERARG, file=sys.stderr)
        sys.exit(1)

    solvernames = [x.strip() for x in sys.argv[pos + 1].split(",")]
    for s in solvernames:
        if s not in allsolvers:
            print("Error: the solver '%s' is not known, please"
                  " use one of: %s" % (s, ", ".join(list(allsolvers.keys()))), file=sys.stderr)
            sys.exit(1)
        solversubsetnames.add(s)
    del sys.argv[pos:pos + 2]

if solversubsetnames:
    # Ensure MipWrapper is included if needed
    if 'CPLEX' in solversubsetnames or 'Gurobi' in solversubsetnames or \
            'SCIP' in solversubsetnames:
        solversubsetnames.add("MipWrapper")

    # Ensure SatWrapper is included if needed
    if 'MiniSat' in solversubsetnames or 'Walksat' in solversubsetnames:
        solversubsetnames.add("SatWrapper")

    extensions = [allsolvers[s] for s in solversubsetnames]


long_desc = """Numberjack is a modelling package written in Python for
constraint programming and combinatorial optimization. Python benefits from a
large and active programming community, Numberjack is therefore a perfect tool
to embed CP technology into larger applications. It is designed to support a
number of efficient underlying C/C++ solvers seamlessly and efficiently.
"""

lic = "License :: OSI Approved :: " \
    "GNU General Public License v2 or later (GPLv2+)"

setup(
    name='Numberjack',
    version='1.2.0',
    author='Numberjack Developers',
    packages=['Numberjack', 'Numberjack.solvers'],
    ext_modules=extensions,
    author_email='numberjack.support@gmail.com',
    url='http://numberjack.ucc.ie/',
    license=lic,
    description='A Python platform for combinatorial optimization.',
    long_description=long_desc,
    cmdclass={'build_ext': njbuild_ext},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Environment :: MacOS X",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: C",
        "Programming Language :: C++",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering",
        lic,
    ],
)
