import os
import shutil
import click
import hashlib
import filecmp


def archive_folders(source_folders, destination_folder,
                    dry=False, verbose=False):
    for folder in source_folders:  # Iterate folders
        if os.path.isdir(folder):
            folder_size = get_file_or_folder_size(folder)
            folder_size = get_human_readable_size(folder_size)
            click.echo("[{}] {}".format(
                folder_size.rjust(10),
                folder))
    if dry:
        click.echo("This is just a dry run. Nothing is executed.")
    else:
        question_text = "Target directory:\n"+destination_folder+"\n\n" \
                        "Do you want to archive (copy to the destination " \
                        "directory and delete the source) all the above " \
                        "folders? "
        user_response = ask_yes_no_question(question_text)
        if user_response:
            for folder in source_folders:
                target_path = os.path.join(destination_folder,
                                           os.path.basename(folder))
                # copy folder
                message = "Copy to " + target_path
                styled_message = click.style(message, fg='bright_green')
                click.echo(styled_message)
                copy_folder(folder, target_path, verbose=True)
                # check all files and directories
                message = "Check all files and directories are identical... "
                styled_message = click.style(message, fg='bright_green')
                click.echo(styled_message)
                compare_result = compare_folders(folder, target_path)
                if compare_result:
                    print("Everything is identical.")
                else:
                    print("Folders are not identical.")
            # delete source
            question_text = "Do you want to delete the source folders now?"
            user_response = ask_yes_no_question(question_text)
            if user_response:
                for folder in source_folders:
                    delete_files_and_folders(folder)


def copy_folder(source_folder, destination_folder, verbose=False):
    try:
        shutil.copytree(source_folder, destination_folder)
        if verbose:
            print("Folder '{}' copied to '{}'.".format(source_folder,
                                                       destination_folder))
    except Exception as e:
        if verbose:
            print(f"Error copying folder: {e}")


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


def get_file_or_folder_size(path):
    if os.path.isfile(path):
        # If it's a file, get its size directly
        size_bytes = os.path.getsize(path)
    elif os.path.isdir(path):
        # If it's a directory, sum up the sizes of all files in it
        size_bytes = sum(os.path.getsize(os.path.join(root, file))
                         for root, dirs, files in os.walk(path)
                         for file in files)
    return size_bytes


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


def calculate_md5(file_path):
    """Calculate the MD5 hash for a file."""
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as file:
        # Read the file in chunks to avoid memory issues with large files
        for chunk in iter(lambda: file.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def compare_folders(folder1, folder2):
    """Compare two folders recursively using MD5 hash values."""
    dcmp = filecmp.dircmp(folder1, folder2)
    identical = True

    # Check for common files and calculate MD5 hash for each
    for common_file in dcmp.common_files:
        file1_path = os.path.join(folder1, common_file)
        file2_path = os.path.join(folder2, common_file)

        hash1 = calculate_md5(file1_path)
        hash2 = calculate_md5(file2_path)

        if hash1 != hash2:
            print(f"Files {file1_path} and {file2_path} are not identical.")
            identical = False

    # Recursively compare subdirectories
    for subdirectory in dcmp.common_dirs:
        subdirectory1 = os.path.join(folder1, subdirectory)
        subdirectory2 = os.path.join(folder2, subdirectory)
        if not compare_folders(subdirectory1, subdirectory2):
            identical = False

    return identical
