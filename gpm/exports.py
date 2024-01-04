import os
import sys
from pathlib import Path
import shutil
import click
import subprocess
import string
import random
from gpm.helper import get_gpmdata_path
import xtarfile as tarfile
from tqdm import tqdm
import hashlib


def check_export_directory(export_folder):
    # Check if the folder exists
    if not os.path.exists(export_folder):
        # If not, create the folder
        os.makedirs(export_folder)
        print(f"Folder '{export_folder}' created.")
    else:
        print(f"Folder '{export_folder}' already exists.")
        sys.exit()


def get_htaccess_path():
    """
    Return the path of htaccess
    """
    config_path = os.path.join(get_gpmdata_path(), "config/htaccess")
    return config_path


def htpasswd_create_user(export_folder, url, username,
                         app=None, raw_export=False):
    """Create the new user in the target directory with password"""
    export_base_path = Path(export_folder).parent.absolute()
    if os.path.exists(os.path.join(export_base_path, ".htpasswd")):
        shutil.copy(os.path.join(export_base_path, ".htpasswd"),
                    os.path.join(export_folder, ".htpasswd"))
        password = generate_password()
        cmd = " ".join(["htpasswd", "-b",
                        os.path.join(export_folder, ".htpasswd"),
                        username, password])
        subprocess.run(cmd, shell=True)
        click.echo()
        click.echo(click.style("Create new user for export directory:",
                               fg='bright_green'))
        click.echo("Directory:\t" + export_folder)
        if app and not raw_export:
            if app in ["RNAseq", "tRNAseq", "mRNAseq", "3mRNAseq"]:
                app = "RNAseq"
            click.echo("".join(["URL:\t", url,
                                "/3_Reports/analysis/Analysis_Report_",
                                app, ".html"]))
        else:
            click.echo("URL:\t" + url)
        click.echo("user:\t" + username)
        click.echo("password:\t" + password)
    else:
        click.echo("Skip setting htpasswd")


def create_user(export_folder, export_URL, username):
    htpasswd_create_user(export_folder, export_URL, username, None)


def generate_password():
    source = string.ascii_letters + string.digits
    result_str = ''.join((random.choice(source) for i in range(12)))
    return result_str


def export_empty_folder(export_URL, export_dir, username):
    # Add htaccess
    # data_dir = os.path.join(os.path.dirname(__file__), "data")
    htaccess_path = get_htaccess_path()
    with open(htaccess_path) as f1:
        contents = [le.rstrip() for le in f1.readlines()]
    for i, line in enumerate(contents):
        if "GPM_TITLE_NAME" in line:
            contents[i] = line.replace("GPM_TITLE_NAME",
                                       os.path.basename(export_dir))

    with open(os.path.join(export_dir, ".htaccess"), "w") as f2:
        for line in contents:
            print(line, file=f2)
    # Create user
    htpasswd_create_user(export_dir,
                         os.path.join(export_URL,
                                      os.path.basename(export_dir)),
                         username.lower(), None)


def relpath(path_file):
    base_dirs = ["/mnt/nextgen", "/mnt/nextgen2", "/mnt/nextgen3"]
    for base_dir in base_dirs:
        if path_file.startswith(base_dir):
            # Getting the relative path of the directory
            rel_path = os.path.relpath(path_file, base_dir)
            print("rel_path: " + rel_path)
            path_file = os.path.join("/", rel_path)
            print("path_file: " + path_file)
    return path_file


def tar_exports(export_folder, dry_run, gzip, same_server=False, prefix=""):
    if export_folder == ".":
        export_folder = os.getcwd()
    export_folder = export_folder.rstrip("/")
    name = os.path.basename(export_folder)
    compressed_folder = os.path.join(export_folder, "compressed_tars")
    # Create compressed_tar folder
    if not os.path.exists(compressed_folder):
        click.echo(click.style("Create the folder:", fg='bright_green'))
        click.echo(compressed_folder)
        if not dry_run:
            os.makedirs(compressed_folder)
    # Tar each folder
    for filename in os.listdir(export_folder):
        if filename.startswith("."):
            continue
        path_file = os.path.join(export_folder, filename)
        if gzip:
            tarfile = os.path.join(compressed_folder,
                                   name+"_" + filename + ".tar.gz")
        else:
            tarfile = os.path.join(compressed_folder,
                                   name+"_" + filename + ".tar")
        # print("path_file: " + path_file)
        if os.path.islink(path_file):
            path_file = os.readlink(path_file)
            # print("path_file link: " + path_file)
            # path_file = relpath(path_file)
            print("path_file link: " + path_file)
        if os.path.isdir(path_file) and filename != "compressed_tars":
            click.echo(click.style("Tar the folder:", fg='bright_green'))
            click.echo(path_file + click.style(" => ",
                                               fg='bright_green') + tarfile)
            if not dry_run:
                if same_server:
                    tar_folder(path_file, tarfile, gzip)
                else:
                    tar_folder_from_project(path_file, tarfile, gzip, prefix=prefix)


def tar_folder(input_folder, output_tar, gzip):
    if gzip:
        tar_mode = 'w:gz'
    else:
        tar_mode = 'w'
    if os.path.exists(output_tar):
        click.echo(output_tar + " exists.")
    else:
        with tarfile.open(output_tar, tar_mode) as tar:
            # Set up the tqdm progress bar
            total_size = 0
            for root, dirs, files in os.walk(input_folder, followlinks=True):
                for name in files:
                    full_path = os.path.join(root, name)
                    total_size += get_size(full_path)
            progress_bar = tqdm(total=total_size,
                                desc='Creating tar archive',
                                unit='B', unit_scale=True)
            for root, dirs, files in os.walk(input_folder, followlinks=True):
                for name in files:
                    full_path = os.path.join(root, name)
                    arcname = os.path.relpath(full_path, input_folder)
                    # Add the file or directory to the tar archive
                    tar.add(full_path, arcname=arcname)
                    # Update the progress bar by the size
                    progress_bar.update(get_size(full_path))
            # Close the progress bar
            progress_bar.close()
        save_md5_to_file(output_tar)


def tar_folder_from_project(input_folder, output_tar, gzip, prefix):
    if gzip:
        tar_mode = 'w:gz'
    else:
        tar_mode = 'w'
    if os.path.exists(output_tar):
        click.echo(output_tar + " exists.")
    else:
        with tarfile.open(output_tar, tar_mode) as tar:
            # Set up the tqdm progress bar
            total_size = 0
            for root, dirs, files in os.walk(input_folder):
                for name in files+dirs:
                    full_path = os.path.join(root, name)
                    if os.path.islink(full_path):
                        print("islink")
                        softlink_path = os.readlink(full_path)
                        print(softlink_path)
                        if softlink_path.startswith(prefix):
                            softlink_path = softlink_path.replace(prefix, "")
                        print(softlink_path)
                        # linked_name = os.basename(full_path)
                        total_size += get_size(softlink_path)
                    else:
                        total_size += get_size(full_path)
            sys.exit()
            progress_bar = tqdm(total=total_size,
                                desc='Creating tar archive',
                                unit='B', unit_scale=True)
            for root, dirs, files in os.walk(input_folder, followlinks=True):
                for name in files:
                    full_path = os.path.join(root, name)
                    arcname = os.path.relpath(full_path, input_folder)
                    # Add the file or directory to the tar archive
                    tar.add(full_path, arcname=arcname)
                    # Update the progress bar by the size
                    progress_bar.update(get_size(full_path))
            # Close the progress bar
            progress_bar.close()
        save_md5_to_file(output_tar)


def tar_folder_from_project2(input_folder, output_tar_gz, project_path):
    def search_name_in_path(query_name):
        res = None
        for root, dirs, files in os.walk(project_path):
            for name in dirs+files:
                if query_name == os.basename(name):
                    res = os.path.join(root, name)
                    break
        return res

    if os.path.exists(output_tar_gz):
        click.echo(output_tar_gz + " exists.")
    else:
        with tarfile.open(output_tar_gz, 'w:gz') as tar:
            # Set up the tqdm progress bar
            total_size = 0
            for root, dirs, files in os.walk(input_folder):
                for name in files:
                    full_path = os.path.join(root, name)
                    print(full_path)
                    if os.islink(full_path):
                        print("islink")
                        linked_name = os.basename(full_path)
                        local_linked = search_name_in_path(linked_name)
                        if local_linked:
                            total_size += get_size(local_linked)
                        else:  # Not found in the project
                            print("Not found in the project: " + full_path)
                            continue
                    else:
                        total_size += get_size(full_path)
            progress_bar = tqdm(total=total_size,
                                desc='Creating tar archive',
                                unit='B', unit_scale=True)
            for root, dirs, files in os.walk(input_folder):
                for name in files:
                    full_path = os.path.join(root, name)
                    arcname = os.path.relpath(full_path, input_folder)
                    # Add the file or directory to the tar archive
                    tar.add(full_path, arcname=arcname)
                    # Update the progress bar by the size
                    progress_bar.update(get_size(full_path))
            # Close the progress bar
            progress_bar.close()
        save_md5_to_file(output_tar_gz)


def get_size(path):
    # If path is a symbolic link, get the size of the target file
    if os.path.islink(path):
        target_path = os.path.realpath(path)
        target_stat = os.stat(target_path)
        return target_stat.st_size
    else:
        # If not a symbolic link, get the size directly
        return os.path.getsize(path)


def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            md5.update(chunk)
    return md5.hexdigest()


def save_md5_to_file(file_path):
    md5_checksum = calculate_md5(file_path)
    with open(file_path+".md5", 'w') as md5_file:
        md5_file.write(md5_checksum)
