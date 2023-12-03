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


def get_gpm_config(section, item):
    """
    Return the config from GPMDATA/gpm.config as a dictionary.
    User-defined config (gpm.config.user) has a higher priority
    than default one (gpm.config).
    """
    gpm_config = path.join(get_gpmdata_path(), "gpm.config.user")
    if not path.exists(gpm_config):
        gpm_config = path.join(get_gpmdata_path(), "gpm.config")

    config = configparser.ConfigParser()
    config.read(gpm_config)
    res = config[section][item]
    if "," in res:
        res = [x.strip() for x in res.split(",")]
    return res
