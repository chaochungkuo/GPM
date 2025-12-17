"""
Minimal setup.py for post-install hooks.

NOTE: This file is kept only to handle the post-install data folder copying.
All package metadata (dependencies, version, etc.) is defined in pyproject.toml.
This file will be removed in a future version once a pure pyproject.toml solution
for post-install hooks is available.
"""

import os
import sys
from os import makedirs, path, walk
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


def get_data_files():
    """
    Collect data files from config, demultiplex, processing, and analysis directories.

    Returns:
        list: List of tuples (install_directory, [list of files]) for setuptools data_files
    """
    data_folders = ["config", "demultiplex", "processing", "analysis"]
    project_root = path.dirname(path.abspath(__file__))
    data_files_list = []

    for folder in data_folders:
        source_path = path.join(project_root, folder)
        if path.exists(source_path) and path.isdir(source_path):
            # Collect all files recursively
            files = []
            for root, dirs, filenames in walk(source_path):
                for filename in filenames:
                    # Get relative path from project root
                    file_path = path.join(root, filename)
                    rel_path = path.relpath(file_path, project_root)
                    files.append(rel_path)

            if files:
                # Install to share/gpm/{folder}/
                install_dir = path.join("share", "gpm", folder)
                data_files_list.append((install_dir, files))

    return data_files_list


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
    data_files=get_data_files(),
    cmdclass={
        "install": PostInstallCommand,
    },
)
