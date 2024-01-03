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


def tar_exports(export_folder, dry_run, same_server=False):
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
        tarfile = os.path.join(compressed_folder, name+"_" + filename + ".tar")
        print("path_file: " + path_file)
        if os.path.islink(path_file):
            path_file = os.readlink(path_file)
            print("path_file link: " + path_file)
            path_file = relpath(path_file)
            print("path_file link: " + path_file)
        if os.path.isdir(path_file) and filename != "compressed_tars":
            click.echo(click.style("Tar the folder:", fg='bright_green'))
            click.echo(path_file + click.style(" => ",
                                               fg='bright_green') + tarfile)
            if not dry_run:
                if same_server:
                    tar_folder(path_file, tarfile)
                else:
                    tar_folder(path_file, tarfile)


def tar_dir(path, tar_name):
    # cmd = " ".join(["tar", "-hcf", tar_name, "-C", os.path.dirname(path),
    #                 "--absolute-names", path])
    # TODO testing tar link dir
    # tar_basepath = tar_name.replace(".tar", "")
    if os.path.exists(tar_name):
        click.echo(tar_name + " exists.")
    else:
        cmd = " ".join(["tar cfh - -C", os.path.dirname(path), path,
                        "-P | pv -s $(du -sb -L ", path,
                        "| awk '{print $1}') >",
                        tar_name])
        # print(cmd)
        # subprocess.run(cmd, shell=True,
        #                stdout=subprocess.PIPE,
        #                stderr=subprocess.PIPE, text=True)
        result = subprocess.run(cmd, shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True)
        if result.stderr:
            print("STDOUT:")
            print(result.stdout)
            print("\nSTDERR:")
            print(result.stderr)

            # for filename in os.listdir(path):
            #     if filename.startswith("."):
            #         continue
            #     path_file = os.path.join(path, filename)
            #     tarfile = os.path.join(path, tar_basepath+"_" + filename + ".tar")
            #     print("path_file: " + path_file)
            #     if os.path.islink(path_file):
            #         path_file = os.readlink(path_file)
            #         print("path_file link: " + path_file)
            #         path_file = relpath(path_file)
            #         print("path_file link: " + path_file)

            #         cmd = " ".join(["tar cfh - -C", os.path.dirname(path_file),
            #                         path_file,
            #                         "-P | pv -s $(du -sb -L ", path_file,
            #                         "| awk '{print $1}') >",
            #                         tarfile])
            #         subprocess.run(cmd, shell=True,
            #                        stdout=subprocess.PIPE,
            #                        stderr=subprocess.PIPE, text=True)
            #         run_md5sum(tar_name)
        # subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # subprocess.call(cmd, shell=True,
        #                                  stdout=subprocess.DEVNULL)
        else:
            run_md5sum(tar_name)
        # subprocess.call(cmd, shell=True,
        #                                  stdout=subprocess.DEVNULL)


def run_md5sum(tar_name):
    cmd = " ".join(["md5sum", tar_name, ">", tar_name+".md5"])
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def tar_dir_distanced(path, tar_name):
    cmd = " ".join(["tar cfh - -C", os.path.dirname(path), path,
                    "-P | pv -s $(du -sb -L ", path, "| awk '{print $1}') >",
                    tar_name])
    print(cmd)
    # subprocess.run(cmd, shell=True,
    #                stdout=subprocess.PIPE,
    #                stderr=subprocess.PIPE, text=True)
    result = subprocess.run(cmd, shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    print("STDOUT:")
    print(result.stdout)
    print("\nSTDERR:")
    print(result.stderr)


def tar_folder(input_folder, output_tar):
    with tarfile.open(output_tar, 'w') as tar:
        # Get all files and directories in the input folder
        files_and_dirs = [os.path.join(input_folder, f)
                          for f in os.listdir(input_folder)]
        # Set up the tqdm progress bar
        total_size = sum(get_size(f) for f in files_and_dirs)
        progress_bar = tqdm(total=total_size,
                            desc='Creating tar archive',
                            unit='B', unit_scale=True)
        for file_or_dir in files_and_dirs:
            # Add the file or directory to the tar archive
            tar.add(file_or_dir,
                    arcname=os.path.relpath(file_or_dir, input_folder))
            progress_bar.update(get_size(file_or_dir))
        # Close the progress bar
        progress_bar.close()


def get_size(path):
    # If path is a symbolic link, get the size of the target file
    if os.path.islink(path):
        target_path = os.path.realpath(path)
        target_stat = os.stat(target_path)
        return target_stat.st_size
    else:
        # If not a symbolic link, get the size directly
        return os.path.getsize(path)
