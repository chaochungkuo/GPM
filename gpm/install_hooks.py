"""
Post-install hooks for GPM package installation.

This module handles the copying of data folders (config, demultiplex, processing, analysis)
to the GPMDATA location during package installation.
"""

import os
import sys
from os import makedirs, path
from shutil import copytree


def get_gpmdata_path():
    """
    Determine the GPMDATA location from environment variable or default to ~/gpmdata.
    
    Returns:
        str: Path to GPMDATA directory
    """
    if os.environ.get("GPMDATA"):
        gpm_data_location = path.expanduser(os.getenv("GPMDATA"))
    else:
        gpm_data_location = path.expanduser(path.join(os.getenv("HOME"), "gpmdata"))
    return gpm_data_location


def copy_gpmdata_folders():
    """
    Copy data folders (config, demultiplex, processing, analysis) to GPMDATA location.
    
    This function:
    1. Determines GPMDATA location from environment variable or defaults to ~/gpmdata
    2. Creates the GPMDATA directory if it doesn't exist
    3. Copies config, demultiplex, processing, and analysis folders
    """
    # Get GPMDATA path
    gpm_data_location = get_gpmdata_path()
    
    # Print location for user feedback
    print("GPMDATA folder: " + gpm_data_location)
    sys.stdout.flush()
    
    # Create Data Path if it doesn't exist
    if not path.exists(gpm_data_location):
        makedirs(gpm_data_location)
    
    # Copy data folders
    data_folders = ["config", "demultiplex", "processing", "analysis"]
    
    # Get the project root directory (where setup.py/pyproject.toml is located)
    # This assumes the package is being installed from source
    project_root = path.dirname(path.dirname(path.abspath(__file__)))
    
    for copy_folder in data_folders:
        source_path = path.join(project_root, copy_folder)
        copy_dest_path = path.join(gpm_data_location, copy_folder)
        
        # Only copy if source exists (handles both source and installed package cases)
        if path.exists(source_path):
            copytree(source_path, copy_dest_path, dirs_exist_ok=True)


if __name__ == "__main__":
    # Allow running as a script for testing
    copy_gpmdata_folders()

