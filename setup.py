#!/usr/bin/env python

from distutils.core import setup

setup(name='magnummsh',
      version='2.0.0',
      description='magnum.msh is the mash converter of magnum.fe',
      author='Claas Abert, Florian Bruckner',
      author_email='claas.abert@tuwien.ac.at, florian.bruckner@tuwien.ac.at',
      url='http://github.com/micromagnetics/magnum.msh',
      packages=['magnummsh'],
      scripts=['scripts/magnum.msh', 'scripts/inp2magnum', 'scripts/magnum2inp']
     )
