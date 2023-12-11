import sys
import io
import os
import re
from setuptools import setup, find_packages
from os import path, makedirs
from shutil import copytree

#############################################################
# Get version
#############################################################


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


current_version = find_version("gpm", "__version__.py")

#############################################################
# GPMDATA path
#############################################################

# if the environment variable is set, use it;
# otherwise use the home directory as a default
if os.environ.get('GPMDATA'):
    gpm_data_location = path.expanduser(os.getenv("GPMDATA"))
else:
    gpm_data_location = path.expanduser(path.join(os.getenv("HOME"),
                                                  "gpmdata"))
print("GPMDATA folder: "+gpm_data_location)
sys.stdout.flush()

# Creating Data Path
if not path.exists(gpm_data_location):
    makedirs(gpm_data_location)

#############################################################
# Copying data
#############################################################

data_folders = ["config", "demultiplex", "processing", "analysis"]
# GPM Configs

for copy_folder in data_folders:
    copy_dest_path = path.join(gpm_data_location, copy_folder)
    copytree(copy_folder, copy_dest_path, dirs_exist_ok=True)
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

short_description = 'GPM (Genomic Project Manager) is a versatile '
'command-line tool designed for managing and automating bioinformatic '
'workflows.'

setup(
    name='gpm',
    version=current_version,
    author='Chao-Chung Kuo',
    author_email='chao-chung.kuo@rwth-aachen.de',
    description=short_description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/chaochungkuo/GPM',
    packages=find_packages(),
    install_requires=[
        'Click', "pandas", "pyyaml"
    ],
    entry_points={
        'console_scripts': [
            'gpm=gpm.main:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
