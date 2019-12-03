#!/usr/bin/env python

from setuptools import setup
import glob
import os


with open("README.md", "r") as f:
    long_description = f.read()

def get_scripts(generic):
    """ Grab all the scripts in the bin directory.  """
    scripts = []
    if os.path.isdir('bin'):
        scripts = [ fname for fname in glob.glob(os.path.join('bin', '*')) ]
    return scripts

setup(name='specklepy',
      version='0.4.0',
      description='',
      long_description=long_description,
      long_description_content_type="text/markdown",
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Data Processing :: Astronomy',
        'Operating System :: OS Independent',
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
      scripts=get_scripts('specklepy/scripts/*py'),
      package_data={},
      include_package_data=True,
      zip_safe=False)
