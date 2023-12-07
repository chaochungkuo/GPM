import sys
from os import path, getenv, environ
import datetime
import click
import configparser
from gpm import CONFIG_LIST


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


def remove_end_slash(path):
    """
    Remove the ending slash in the given path.

    :param path: A path.
    :type path: str
    :return: the modified path
    :rtype: str
    """
    if path.endswith("/"):
        path = path.rstrip("/")
    return path


#################################################################
# Config relevant
#################################################################
def get_config_values(config_name):
    """
    Return a dictionary of the key value pairs in the defined config file.
    User-defined config (*.config.user) has a higher priority
    than default one (*.config).
    """
    config_path = path.join(get_gpmdata_path(), "config/"+config_name+".user")
    if not path.exists(config_path):
        config_path = path.join(get_gpmdata_path(), "config/"+config_name)
    config = configparser.ConfigParser()
    config.read(config_path)
    combined_dict = {}
    for section_name in config.sections():
        section_dict = dict(config.items(section_name))
        combined_dict.update(section_dict)
    return combined_dict


def get_config_value(config_name, section, item):
    """
    Return value of the defined item in the defined config file.
    User-defined config (*.config.user) has a higher priority
    than default one (*.config).
    """
    config_path = path.join(get_gpmdata_path(), "config/"+config_name+".user")
    if not path.exists(config_path):
        config_path = path.join(get_gpmdata_path(), "config/"+config_name)
    config = configparser.ConfigParser()
    config.read(config_path)
    res = config[section][item]
    if "," in res:
        res = [x.strip() for x in res.split(",")]
    return res


def get_gpm_config(section, item):
    """
    Return the config from GPMDATA/gpm.config as a dictionary.
    User-defined config (gpm.config.user) has a higher priority
    than default one (gpm.config).
    """
    config_name = "gpm.config"
    res = get_config_value(config_name, section, item)
    return res


def get_environment_config(section, item):
    """
    Return the config from GPMDATA/environment.config as a dictionary.
    User-defined config (environment.config.user) has a higher priority
    than default one (environment.config).
    """
    config_name = "environment.config"
    res = get_config_value(config_name, section, item)
    return res


def get_dict_from_configs():
    combined_dict = {}
    for sel_config in CONFIG_LIST:
        config_dict = get_config_values(sel_config)
        combined_dict.update(config_dict)
    return combined_dict


def replace_variables_by_dict(line, input_dict):
    for key, value in input_dict.items():
        if key.upper() in line:
            line = line.replace(key.upper(), value)
    return line


def check_analysis_name(analysis_dict, analysis_name):
    all_names = []
    for k, g in analysis_dict.items():
        all_names = all_names + list(g.keys())
    if analysis_name not in all_names:
        click.echo("Please choose an analysis from the list below")
        click.echo(all_names)
        sys.exit()
