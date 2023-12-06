from os import path
import configparser
from gpm import CONFIG_LIST
from gpm.helper import get_gpmdata_path


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
