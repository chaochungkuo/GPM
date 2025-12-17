"""
Minimal setup.py for post-install hooks.

NOTE: This file is kept only to handle the post-install data folder copying.
All package metadata (dependencies, version, etc.) is defined in pyproject.toml.
This file will be removed in a future version once a pure pyproject.toml solution
for post-install hooks is available.
"""

import os
import sys
from os import makedirs, path
from shutil import copytree

from setuptools import setup
from setuptools.command.install import install


def get_gpmdata_path():
    """Determine the GPMDATA location from environment variable or default to ~/gpmdata."""
    if os.environ.get("GPMDATA"):
        gpm_data_location = path.expanduser(os.getenv("GPMDATA"))
    else:
        gpm_data_location = path.expanduser(path.join(os.getenv("HOME"), "gpmdata"))
    return gpm_data_location


def copy_gpmdata_folders():
    """Copy data folders to GPMDATA location."""
    gpm_data_location = get_gpmdata_path()
    print("GPMDATA folder: " + gpm_data_location)
    sys.stdout.flush()

    if not path.exists(gpm_data_location):
        makedirs(gpm_data_location)

    data_folders = ["config", "demultiplex", "processing", "analysis"]
    # Get project root (where this setup.py file is located)
    project_root = path.dirname(path.abspath(__file__))

    for copy_folder in data_folders:
        source_path = path.join(project_root, copy_folder)
        copy_dest_path = path.join(gpm_data_location, copy_folder)
        if path.exists(source_path) and path.isdir(source_path):
            copytree(source_path, copy_dest_path, dirs_exist_ok=True)
        else:
            print(f"Warning: Source folder '{source_path}' not found. Skipping.")
            sys.stdout.flush()


class PostInstallCommand(install):
    """Custom install command to run post-install hooks."""

    def run(self):
        """Run the standard install, then copy GPMDATA folders."""
        install.run(self)
        # copy_gpmdata_folders()


setup(
    cmdclass={
        "install": PostInstallCommand,
    },
)
