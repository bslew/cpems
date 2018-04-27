#!/usr/bin/env python

from distutils.core import setup

setup(name='pyCPEDScommonFunctions',
      version='1.0',
      description='A bundle of common routines from cpems package',
      author='Bartosz Lew',
      author_email='blew@astro.umk.pl',
      url='',
      package_dir = {'': ''},
      packages = ['pyCPEDScommonFunctions'],
      scripts=['pyCPEDScommonFunctions/confidenceRange.py','pyCPEDScommonFunctions/cal2jd.py',]
     )

#       py_modules=['RadiometerData.RPG_tau','RadiometerData.RPG_Tatm'],


