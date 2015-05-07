from setuptools import setup
from setuptools import find_packages
import os
import jacquard

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as filename:
        return filename.read()

setup(name='jacquard',
      version=jacquard.__version__,
      description=('Command-line tools to expedite analysis of '
                   'Variant Call Format (VCF) files.'),
      long_description=(read('README.rst') + '\n\n' +
                        read('CHANGELOG.rst') + '\n\n' +
                        read('AUTHORS.rst')),
      url='https://github.com/umich-brcf-bioinf/Jacquard',
      author='University of Michigan Bioinformatics Core',
      author_email='bfx-jacquard@umich.edu',
      license='Apache',
      packages=find_packages(exclude=['test*']),
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: Apache Software License',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'],
      keywords='VCF bioinformatic exome-seq DNA-seq variant-call-format',
      install_requires=['natsort'],
      entry_points={'console_scripts': ['jacquard=jacquard.jacquard:main']},
      test_suite='nose.collector',
      tests_require=['nose', 'testfixtures', 'numpy'],
      zip_safe=False)
