from os import path, getenv, environ
import configparser


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


def load_gpm_config():
    """
    Return the config from GPMDATA/gpm.config as a dictionary
    """
    gpmdata_path = get_gpmdata_path()
    Config = configparser.ConfigParser()
    Config.read(path.join(gpmdata_path, "gpm.config"))
    return Config
