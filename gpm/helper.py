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
    # Retrieve values from the configuration file
    section = config['Section']
    key1 = section.get('key1')
    key2 = section.getint('key2')  # Convert to int
    key3 = section.getboolean('key3')  # Convert to boolean
    
    # Parse list from a comma-separated string
    list_key = section.get('list_key').split(', ')
    
    # Return a dictionary with the loaded values
    config_values = {
        'key1': key1,
        'key2': key2,
        'key3': key3,
        'list_key': list_key,
    }
    
    return config_values
    return Config
