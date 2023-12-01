from setuptools import setup, find_packages
# import os
from gpm.version import version
from glob import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

short_description = 'The Genomic Project Manager is a powerful tool designed '
'to streamline and automate bioinformatic workflows within a '
'core facility.'

setup(
    name='gpm',
    version=version,
    author='Chao-Chung Kuo',
    author_email='chao-chung.kuo@rwth-aachen.de',
    description=short_description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/chaochungkuo/GPM',
    packages=find_packages(),
    install_requires=[
        'Click', 'pandas'
    ],
    entry_points={
        'console_scripts': [
            'gpm=gpm.main:main',
        ],
    },
    data_files=[
        ('gpm/config', glob('config/*')),
        ('gpm/analysis', glob('analysis/*')),
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
