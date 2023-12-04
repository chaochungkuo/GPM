#!/usr/bin/env python
import sys
from os import path, makedirs, listdir
from shutil import copyfile
from helper import get_gpmdata_path


def prebuild_script():
    gpm_data_location = get_gpmdata_path()
    print("GPMDATA folder: "+gpm_data_location)
    create_data_path(gpm_data_location)

    installation_path = path.dirname(__file__)
    print(installation_path)
    sys.stdout.flush()
    # if path.exists('gpm/demultiplex'):
    #     rmtree('gpm/demultiplex')


def create_data_path(gpm_data_location):
    # Creating Data Path
    if not path.exists(gpm_data_location):
        makedirs(gpm_data_location)
    # GPM Configs
    config_dir = path.join(path.dirname(path.dirname(__file__)), "config")
    for file in listdir(config_dir):
        fn = path.basename(file)
        copyfile(path.join(config_dir, fn),
                 path.join(gpm_data_location, fn))
        # User defined Configs
        # userconfig = open(path.join(gpm_data_location,fn+".user"), "w")
        # with open(path.join(config_dir,fn)) as f1:
        #     for line in f1.readlines():
        #         print("# "+line, file=userconfig, end="")
        # userconfig.close()


if __name__ == "__main__":
    prebuild_script()
