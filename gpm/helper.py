import sys
from os import path, getenv, environ
import datetime
import click
from gpm.configs import get_gpm_config


def get_gpmdata_path():
    """
    Get the GPMDATA path from the environment.
    """
    if environ.get('GPMDATA'):
        gpm_data_location = path.expanduser(getenv("GPMDATA"))
    else:
        gpm_data_location = path.expanduser(path.join(getenv("HOME"),
                                                      "gpmdata"))
    return gpm_data_location


def check_project_name(name):
    """
    Check the name of the project with the pattern
    YYMMDD_Name1_Name2_Institute_Application
    """
    split_name = name.split("_")
    if len(split_name) != 5:
        print("Error: Please follow the pattern below:")
        print("YYMMDD_Name1_Name2_Institute_Application")
        sys.exit()
    # date
    sequencing_date = split_name[0]
    try:
        datetime.datetime.strptime(sequencing_date, "%y%m%d")
    except ValueError:
        click.echo("Incorrect date string format. It should be YYMMDD")
        sys.exit()
    # app
    app = split_name[4]
    if app not in get_gpm_config("GPM", "APPLICATIONS"):
        click.echo("Unsupported application. Please take one from below:")
        click.echo(", ".join(get_gpm_config("GPM", "APPLICATIONS")))
        sys.exit()
