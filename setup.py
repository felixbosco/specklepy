#!/usr/bin/env python

from setuptools import setup, find_namespace_packages
import glob
import os


def long_description():
    with open("README.md", "r") as f:
        long_description = f.read()
    return long_description

def find_scripts():
    """Grab all the scripts in the bin directory."""
    scripts = []
    if os.path.isdir('bin'):
        scripts = glob.glob(os.path.join('bin', '*'))
    return scripts

def find_packages():
    return find_namespace_packages(exclude=['*build/*', '*data/*', '*deprecated/*'])



setup(name='specklepy',
      version='0.5.1',
      description='Specklepy Holographic Data Reduction',
      long_description=long_description(),
      long_description_content_type="text/markdown",
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'Intended Audience :: Science/Research',
        'Topic :: Data Processing :: Astronomy',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      keywords='',
      url='https://github.com/felixbosco/specklepy',
      author='Felix Bosco',
      author_email='bosco@mpia.de',
      license='MIT',
      install_requires=[
          'astropy',
          'photutils',
      ],
      packages=find_packages(),
      scripts=find_scripts(),
      # package_data={},
      include_package_data=True,
      zip_safe=False)
