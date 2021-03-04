#!/usr/bin/env python

from distutils.core import setup, Extension
import os

os.environ["CC"] = "c++" 
os.environ["CXX"] = "c++"

# # Common flags for both release and debug builds.
# extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
# extra_compile_args += ["-std=c++11", "-Wall", "-Wextra"]
# if _DEBUG:
#     extra_compile_args += ["-g3", "-O0", "-DDEBUG=%s" % _DEBUG_LEVEL, "-UNDEBUG"]
# else:
#     extra_compile_args += ["-DNDEBUG", "-O3"]
    
    
#     libraries=['CPEDS','Mscscore','Mscsfn','hdf5'],
cpedsRotation = Extension(
    'pyCPEDScommonFunctions/cpedsRotation',
    sources=['pyCPEDScommonFunctions/cpedsRotation.cpp'],
    include_dirs=['/usr/local/include/cpems', 
                '/usr/lib64/python3.6/site-packages/numpy/core/include/numpy',
                  '/usr/lib/python3.6/dist-packages/numpy/core/include/numpy/',
                  ],
    library_dirs=['/usr/local/lib/cpems'],
    libraries=['nova',  'gsl', 'gslcblas', 'm', 'proj', 
               'QtCore', 'fftw3', 'fftw3l', 'hdf5', 'CGAL', 'gmp','cfitsio', 'CPEDS', 'Mscsfn', 
               'Mscscore', 'Mscsplot', 'MscsWMAP', 'armadillo', 
               'gsl', 'gslcblas', 'm', 'proj', 'QtCore', 'fftw3', 'ccSHT3', 'novas', 'velKB', 'slaRefr', 
               'fftw3l', 'hdf5', 'CGAL', 'gmp', 'cfitsio', 'cpgplot', 'armadillo'],
    language='C++',
    )

setup(name='pyCPEDScommonFunctions',
      version='1.0',
      description='A bundle of common routines from cpems package',
      author='Bartosz Lew',
      author_email='blew@astro.umk.pl',
      url='',
      package_dir = {'': ''},
      packages = ['pyCPEDScommonFunctions'],
      scripts=['pyCPEDScommonFunctions/confidenceRange.py',
               'pyCPEDScommonFunctions/cal2jd.py',
               'pyCPEDScommonFunctions/jd2cal.py', 
               'pyCPEDScommonFunctions/join_interpolate.py'],
      ext_modules=[cpedsRotation]
     )

#       py_modules = ['pyCPEDScommonFunctions.libcpedsRotation.so'],
#       py_modules=['RadiometerData.RPG_tau','RadiometerData.RPG_Tatm'],


