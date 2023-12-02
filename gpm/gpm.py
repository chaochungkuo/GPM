"""
GPM - Python library for managing genomics projects.
"""
from collection import OrderedDict
import configparser
from gpm.helper import get_gpmdata_path

tags_GPM = OrderedDict([("Project", ["date", "name1", "name2", "institute",
                                     "application"]),
                        ("Raw data", ["bcl_path"]),
                        ("Demultiplexing", ["fastq_path", "fastq_qc_path",
                                            "demultiplex_method"]),
                        ("Processing", ["project_path",
                                        "processing_method",
                                        "processing_qc_path"]),
                        ("Analyses", ["analysis_path",
                                      "analysis_types"]),
                        ("Export", ["export_URL",
                                    "export_user",
                                    "export_password"]),
                        ("History", ["logs"])])


class GPM():
    """A class for a genomic project."""
    def __init__(self):
        """
        Initiate a project.

        :return: None
        """
        self.profile = OrderedDict()
        for section in tags_GPM.keys():
            self.profile[section] = OrderedDict()
            for tag in tags_GPM[section]:
                self.profile[section][tag] = ""

    def load_project_config_file(self, filepath):
        """
        Read a project config file (project.ini)

        :param filepath: Path to read the project config.
        :type filepath: str
        :return: None
        """
        config = configparser.ConfigParser()
        config.read(filepath)
        # Retrieve values from the configuration file
        for section in tags_GPM.keys():
            section_dict = config[section]
            for tag in tags_GPM[section]:
                value = section_dict.get(tag)
                if "," in value:
                    value = [x.strip() for x in value.split(",")]
                elif value == "":
                    continue
                self.profile[section][tag] = value

    def write_project_config_file(self, filepath):
        """
        Write a project config file (project.ini)

        :param filepath: Path to write the project config.
        :type filepath: str
        :return: None
        """
        # Create a configparser object
        config = configparser.ConfigParser()
        # Iterate through the dictionary and add sections and key-value pairs
        for section, options in self.profile.items():
            config.add_section(section)
            for key, value in options.items():
                config.set(section, key, str(value))

        # Write the configuration to the file
        with open(filepath, 'w') as config_file:
            config.write(config_file)

    def copy_file(self, source, target):
        """
        Copy a file from source to target with modifying certain key words

        :param source: Path of the source file.
        :type source: str
        :param target: Path of the target file.
        :type target: str
        :return: None
        """
        with open(source, 'r') as input_file, open(target, 'w') as output_file:
            for line in input_file:
                for section, options in self.profile.items():
                    for tag, value in options.items():
                        if tag.upper() in line:
                            line = line.replace(tag.upper(), value)
                output_file.write(line)

    def load_gpm_config(self):
        """
        Read the GPM config file (gpm.config or gpm.config.user)

        :return: None
        """
        gpmdata = get_gpmdata_path()
        