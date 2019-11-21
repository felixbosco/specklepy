from setuptools import setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(name='HoloPy',
      version='0.0.3',
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
      url='https://github.com/felixbosco/holopy',
      author='Felix Bosco',
      author_email='bosco@mpia.de',
      license='MIT',
      install_requires=[
          'astropy',
          'photutils',
      ],
      package_data={},
      include_package_data=True,
      zip_safe=False)
