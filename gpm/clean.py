import fnmatch
# import glob
import os
import shutil
# import re
import click
from gpm.helper import get_gpm_config


def clean_folders(target_folders, show_each_file, dry=False):
    regex_patterns = get_gpm_config("CLEAN", "PATTERNS")
    paths_to_be_cleaned = []
    for folder in target_folders:  # Iterate folders
        if os.path.isdir(folder):
            folder_size = get_file_or_folder_size(folder, show_each_file)
            if not folder_size:
                continue
            folder_size = get_human_readable_size(folder_size)
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
                total_size = get_human_readable_size(total_size)
                click.echo(click.style("[{}/{}] {}".format(
                    total_size.rjust(10),
                    folder_size.rjust(10),
                    folder), fg='bright_green'))
                if show_each_file:
                    for matching_file in matching_files:
                        readable_size = get_human_readable_size(
                            each_size[matching_file]
                            )
                        click.echo("[{}] {}".format(readable_size.rjust(10),
                                                    matching_file))
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
            for file_or_folder in fnmatch.filter(files + dirs, p):
                matching_files.append(os.path.join(root, file_or_folder))
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
