"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
import os

here = os.path.abspath(os.path.dirname(__file__))
bindir = os.path.join(here, "bin/")

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
	long_description = f.read()

setup(name='gpseqc',
	version='2.3.5',
	description='GPSeq-based centrality estimation.',
	long_description=long_description,
	long_description_content_type='text/markdown',
	url='https://github.com/ggirelli/gpseqc',
	author='Gabriele Girelli',
	author_email='gabriele.girelli@scilifelab.se',
	license='MIT',
	classifiers=[
		'Development Status :: 4 - Beta',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 3 :: Only',
	],
	keywords='DNA biology cell centrality nucleus genome region bed',
	packages=find_packages(),
	install_requires=['ggc>=0.0.3', 'joblib>=0.11', 'matplotlib>=2.2.2',
		'numpy>=1.14.2', 'pandas>=0.22.0', 'pybedtools>=0.7.10', 'pyemd>=0.5.1',
		'scipy>=1.0.1'],
	scripts=['bin/gpseqc_estimate', 'bin/gpseqc_compare'],
	test_suite="nose.collector",
	tests_require=["nose"],
)
