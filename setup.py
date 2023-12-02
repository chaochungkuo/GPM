from setuptools import setup, find_packages
from os import path, getenv, makedirs, listdir, environ
from gpm.version import version
from glob import glob
import shutil

#############################################################
# Creating Data Path
#############################################################

# if the environment variable is set, use it;
# otherwise use the home directory as a default
if environ.get('GPMDATA'):
    gpm_data_location = path.expanduser(getenv("GPMDATA"))
else:
    gpm_data_location = path.expanduser(path.join(getenv("HOME"), "gpmdata"))
print("GPMDATA folder: "+gpm_data_location)

# Creating Data Path
if not path.exists(gpm_data_location):
    makedirs(gpm_data_location)
# GPM Configs
config_dir = path.join(path.dirname(__file__), "config")
for file in listdir(config_dir):
    fn = path.basename(file)
    shutil.copyfile(path.join(config_dir, fn),
                    path.join(gpm_data_location, fn))
    # User defined Configs
    # userconfig = open(path.join(gpm_data_location,fn+".user"), "w")
    # with open(path.join(config_dir,fn)) as f1:
    #     for line in f1.readlines():
    #         print("# "+line, file=userconfig, end="")
    # userconfig.close()

#############################################################
# Setup function
#############################################################

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
