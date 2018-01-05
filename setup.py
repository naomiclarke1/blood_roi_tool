#!/usr/bin/env python

import os
from setuptools import setup

setup(name='blood_tools',
      version='0.0.1',
      description='Blood Analysis GUI',
      url='',
      maintainer='SP',
      maintainer_email='sp',
      license='BSD',
      keywords=[],
      packages=[],
      install_requires=[open('requirements.txt').read().strip().split('\n')],
      long_description=(open('README').read() if os.path.exists('README')
                        else ''),
      entry_points='''
        [gui_scripts]
        blood_tools=blood_tools.cli.main:main
      ''',
      zip_safe=False)
