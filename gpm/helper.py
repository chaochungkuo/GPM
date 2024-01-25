import sys
from os import path, getenv, environ
import datetime
import click
import configparser
from gpm import CONFIG_LIST
import pandas as pd


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
    if app not in get_gpm_config("GPM", "GPM_APPLICATIONS"):
        click.echo("Unsupported application. Please take one from below:")
        click.echo(", ".join(get_gpm_config("GPM", "GPM_APPLICATIONS")))
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
    User-defined config (*.ini.user) has a higher priority
    than default one (*.ini).
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
    User-defined config (*.ini.user) has a higher priority
    than default one (*.ini).
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
    Return the config from GPMDATA/gpm.ini as a dictionary.
    User-defined config (gpm.ini.user) has a higher priority
    than default one (gpm.ini).
    """
    config_name = "gpm.ini"
    res = get_config_value(config_name, section, item)
    return res


def get_environment_config(section, item):
    """
    Return the config from GPMDATA/environment.ini as a dictionary.
    User-defined config (environment.ini.user) has a higher priority
    than default one (environment.ini).
    """
    config_name = "environment.ini"
    res = get_config_value(config_name, section, item)
    return res


def get_dict_from_configs():
    combined_dict = {}
    for sel_config in CONFIG_LIST:
        config_dict = get_config_values(sel_config)
        combined_dict.update(config_dict)
    # print(combined_dict)
    return combined_dict


def replace_variables_by_dict(line, input_dict):
    for key, value in input_dict.items():
        if key.upper() in line:
            if key == "rmd_authors":
                authors = [author.strip() for author
                           in value.split("\n")]
                authors = ["  - "+author.strip('"') for author in authors]
                line = "\n".join(authors)+"\n"
            elif key == "rmd_institute_logo":
                # print(line)
                # print(value)
                line = line.replace(key.upper(), value)
                # print(line)
            else:
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


def copy_samplesheet(source_samplesheet, target_samplesheet):
    try:
        df = pd.read_csv(source_samplesheet)
        df.drop(columns=["fastq1", "fastq2", "strandness"], inplace=True)
        df[['label1', 'label2']] = df['sample'].str.split('_', expand=True)
        df.to_csv(target_samplesheet, index=False)
    except Exception:
        warning = "No samplesheet from processing path is available.\n"\
                  "Please generate the samplesheet.csv manually."
        click.echo(click.style(warning,
                   fg='red'))


def append_file_to_another(file1, file2):
    with open(file1, 'r') as source_file:
        source_content = source_file.read()
    # Append the content to the destination file
    with open(file2, 'a') as destination_file:
        destination_file.write(source_content)
