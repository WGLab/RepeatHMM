#!/usr/bin/env python

import os,sys
#if sys.version_info[:2] < (2, 7):
#    raise Exception('This version of gensim needs Python 2.7 or later. ')
#elif sys.version_info[:2] > (3, 0):
#    raise Exception('This version of gensim might have errors in Python 3 or later. ')

import setuptools
from distutils.core import Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

extension_mod = Extension("_UnsymmetricPairAlignment", sources=["bin/RepeatHMM_scripts/UnsymmetricPairAlignment/UnsymmetricPairAlignment_wrap.c", "bin/RepeatHMM_scripts/UnsymmetricPairAlignment/UnsymmetricPairAlignment.c"], swig_opts=['-modern', '-I../include'])

setuptools.setup(
    name="RepeatHMM", # Replace with your own username
    version="2.0",
    author="Qian Liu",
    author_email="liuqianhn@gmail.com",
    description="A computational tool for repeat count estimation from long reads",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/WGLab/RepeatHMM",
    #packages=setuptools.find_packages(),
    #packages=['RepeatHMM_scripts', 'RepeatHMM_scripts.UnsymmetricPairAlignment'],
    #package_dir={'RepeatHMM_scripts': 'bin/RepeatHMM_scripts', "RepeatHMM_scripts.UnsymmetricPairAlignment":'bin/RepeatHMM_scripts/UnsymmetricPairAlignment' },
    #package_data={'RepeatHMM_scripts.UnsymmetricPairAlignment': ['bin/RepeatHMM_scripts/UnsymmetricPairAlignment/_UnsymmetricPairAlignment.so']},
    packages=['RepeatHMM_scripts', 'UnsymmetricPairAlignment'],
    package_dir={'RepeatHMM_scripts': 'bin/RepeatHMM_scripts', "UnsymmetricPairAlignment":'bin/RepeatHMM_scripts/UnsymmetricPairAlignment' },    
    package_data={'UnsymmetricPairAlignment': ['bin/RepeatHMM_scripts/UnsymmetricPairAlignment/_UnsymmetricPairAlignment.so']},
    scripts=['bin/repeatHMM.py'],
    ext_modules=[extension_mod], py_modules=['UnsymmetricPairAlignment'],
    classifiers=[
        "Programming Language :: Python",
        'Intended Audience :: Science/Research',
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: OS Independent",
    ],
)
