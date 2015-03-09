from setuptools import setup
import jacquard

def readme():
    with open('README.md') as readme_file:
        return readme_file.read()

setup(name='jacquard',
      version=jacquard.__version__,
      description=('Command-line tools to expedite analysis of '
                   'Variant Call Format (VCF) files.'),
      long_description=readme(),
      url='https://github.com/umich-brcf-bioinf/Jacquard',
      author='University of Michigan Bioinformatics Core',
      author_email='bfx-jacquard@umich.edu',
      license='Apache',
      packages=['jacquard', 'jacquard.variant_callers'],
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: Apache Software License',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'],
      keywords='VCF bioinformatic exome-seq DNA-seq variant-call-format',
      install_requires=['natsort'],
      entry_points={'console_scripts': ['jacquard=jacquard.jacquard:main']},
      test_suite='nose.collector',
      tests_require=['nose', 'numpy', 'testfixtures'],
      zip_safe=False)
