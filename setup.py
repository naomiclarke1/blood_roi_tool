#!/usr/bin/env python

import os
from setuptools import setup

setup(name='blood_tools',
      version='0.0.1',
      description='Blood Analysis GUI',
      url='',
      maintainer='Sharon Portnoy',
      maintainer_email='shportnoy@gmail.com',
      license='BSD',
      keywords=[],
      packages=['blood_tools'],
      install_requires=[open('requirements.txt').read().strip().split('\n')],
      long_description=(open('README.md').read() if os.path.exists('README.md')
                        else ''),
      entry_points='''
        [gui_scripts]
        blood_tools=blood_tools.cli.main:main
      ''',
      zip_safe=False)
