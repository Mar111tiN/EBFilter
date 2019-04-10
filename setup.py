#!/usr/bin/env python

from distutils.core import setup

setup(name='ebfilter',
      version='1.0.0',
      description='Python tools for filtering somatic mutations using beta-binomial sequencing error model.',
      author='Martin Szyska modified from original EBFilter by Yuichi Shiraishi',
      author_email='martinszyska37@gmail.com',
      url='https://github.com/Mar111tiN/EBFilter/',
      package_dir = {'': 'code'},
      packages=['ebfilter'],
      scripts=['EBrun'],
      license='GPL-3'
     )


