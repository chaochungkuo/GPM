import fnmatch
# import glob
import os
import shutil
from datetime import datetime
# import re
import click
from gpm.helper import get_gpm_config

# Function to clean the path, extract date, and compare
def folder_before_date(file_path, target_date_str):
    try:
        # Extract the base name from the path
        base_name = os.path.basename(file_path)
        # Split the base name by underscores and extract the first part
        date_str = base_name.split('_')[0]
        # Parse the date string to a datetime object
        date_obj = datetime.strptime(date_str, "%y%m%d")
        target_date_obj = datetime.strptime(target_date_str, "%y%m%d")
        # Compare the dates
        return date_obj < target_date_obj
    except (ValueError, IndexError):
        # Return False if the format is invalid or extraction fails
        return False

def clean_folders(target_folders, show_each_file, keep_files, before="", dry=False):
    regex_patterns = get_gpm_config("CLEAN", "PATTERNS")
    paths_to_be_cleaned = []
    total_size_to_be_cleaned = 0
    for folder in target_folders:  # Iterate folders
        if os.path.isdir(folder):
            if os.path.exists(os.path.join(folder, ".keep")):
                # Skip folders with .keep file
                continue
            elif before!="" and not folder_before_date(folder, before):
                # Skip folders which is not before the target date
                continue
            else:
                # Get total size
                folder_size = get_file_or_folder_size(folder, show_each_file)
                if not folder_size:
                    continue
                folder_size_str = get_human_readable_size(folder_size)
                # Find matching files
                matching_files = search_files_by_patterns(folder,
                                                        regex_patterns)
                if matching_files:
                    paths_to_be_cleaned += matching_files
                    total_size = 0
                    each_size = {}

                    for matching_file in matching_files:
                        size_bytes = get_file_or_folder_size(matching_file,
                                                            show_each_file)
                        total_size += size_bytes
                        each_size[matching_file] = size_bytes
                    total_size_to_be_cleaned += total_size
                    total_size_str = get_human_readable_size(total_size)
                    percentage = total_size/folder_size * 100
                    click.echo(click.style("[{}/{}] {} ".format(
                                           total_size_str.rjust(10),
                                           folder_size_str.rjust(10),
                                           f"{percentage:.1f}%".rjust(6)), 
                                           fg='bright_green')+folder)
                    if show_each_file:
                        for matching_file in matching_files:
                            readable_size = get_human_readable_size(
                                each_size[matching_file]
                                )
                            click.echo("[{}] {}".format(readable_size.rjust(10),
                                                        matching_file))
    total = get_human_readable_size(total_size_to_be_cleaned)
    click.echo(click.style("[{}] can be cleaned. ".format(total),
                           bold=True, fg='bright_green'))
    if not paths_to_be_cleaned:
        click.echo("No files/folders match the defined patterns.")
        click.echo(get_gpm_config("CLEAN", "PATTERNS"))
    elif not dry:
        question_text = "Do you want to clean all the above files/folders?" \
                        "(Not revertible)"
        user_response = ask_yes_no_question(question_text)
        if user_response:
            for paths_to_delete in paths_to_be_cleaned:
                delete_files_and_folders(paths_to_delete)


def search_files_by_patterns(root_path, patterns):
    matching_files = []
    for root, dirs, files in os.walk(root_path):
        for p in patterns:
            # matching pattern for basename
            for file_or_folder in fnmatch.filter(files + dirs, p):
                matching_files.append(os.path.join(root, file_or_folder))
            # matching pattern for the absolute path, e.g. nfcore*/results*
            if "/" in p:
                p1 = p.split("/")[0]
                p2 = p.split("/")[1]
                root_base = os.path.basename(root)
                if fnmatch.fnmatch(root_base, p1):
                    for folder in fnmatch.filter(dirs, p2):
                        matching_files.append(os.path.join(root, folder))
                
    matching_files = remove_subpaths(matching_files)
    return matching_files

def get_total_size(files):
    total_size = sum(os.path.getsize(file) for file in files)
    return total_size


def get_human_readable_size(size_bytes):
    # Define the units and their respective sizes
    units = ['B', 'KB', 'MB', 'GB', 'TB']
    size = size_bytes

    # Find the appropriate unit
    for unit in units:
        if size < 1024.0:
            break
        size /= 1024.0

    return "{:.2f} {}".format(size, unit.ljust(2))


def get_file_or_folder_size(path, verbose=False):
    try:
        if os.path.isfile(path):
            # If it's a file, get its size directly
            size_bytes = os.path.getsize(path)
        elif os.path.isdir(path):
            # If it's a directory, sum up the sizes of all files in it
            size_bytes = sum(os.path.getsize(os.path.join(root, file))
                             for root, dirs, files in os.walk(path)
                             for file in files)
        return size_bytes
    except FileNotFoundError as e:
        if verbose:
            click.echo(e)
        return None


def ask_yes_no_question(question):
    while True:
        response = input(question + " (yes/no): ").lower()
        if response == "yes":
            return True
        elif response == "no":
            return False
        else:
            print("Please enter 'yes' or 'no'.")


def delete_files_and_folders(path):
    if os.path.isfile(path):
        try:
            os.remove(path)
            print(f"File deleted: {path}")
        except Exception as e:
            print(f"Error deleting file {path}: {e}")
    elif os.path.isdir(path):
        try:
            shutil.rmtree(path)
            print(f"Folder deleted: {path}")
        except Exception as e:
            print(f"Error deleting folder {path}: {e}")
    else:
        print(f"Path not found: {path}")


def remove_subpaths(paths):
    # Normalize and sort paths to ensure parent paths come before subpaths
    paths = sorted(paths, key=lambda p: os.path.abspath(p))
    result = []

    for path1 in paths:
        is_subpath = False
        for path2 in result:
            # Check if path1 is a subpath (or file) under path2
            if os.path.commonpath(
                [os.path.abspath(path1),
                 os.path.abspath(path2)]) == os.path.abspath(path2):
                is_subpath = True
                break
        if not is_subpath:
            result.append(path1)
    
    return result
