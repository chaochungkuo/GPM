import os
import sys
import click
from collections import OrderedDict
import configparser
from datetime import datetime
from gpm.helper import remove_end_slash, get_gpmdata_path, \
                       check_project_name
from gpm.configs import get_dict_from_configs, replace_variables_by_dict
from gpm import PROJECT_INI_FILE

tags_GPM = OrderedDict([("Project", ["date", "name1", "name2", "institute",
                                     "application", "project.ini",
                                     "project_path"]),
                        ("Raw data", ["bcl_path"]),
                        ("Demultiplexing", ["fastq_path", "fastq_qc_path",
                                            "demultiplex_method"]),
                        ("Processing", ["processing_path",
                                        "processing_method",
                                        "processing_qc_path",
                                        "organism",
                                        "genome_assembly"]),
                        ("Analyses", ["analysis_path",
                                      "analysis_types"]),
                        ("Export", ["export_URL",
                                    "export_user",
                                    "export_password"])])


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
        self.logs = []

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
            if section == "Logs":
                for option in config.options(section):
                    if not config.has_option(section, option):
                        value = config.get(section, option)
                        self.logs.append(value)
            else:
                section_dict = config[section]
                for tag in tags_GPM[section]:
                    value = section_dict.get(tag)
                    if "," in value:
                        value = [x.strip() for x in value.split(",")]
                    elif value == "":
                        continue
                    self.profile[section][tag] = value

    def write_project_config_file(self):
        """
        Write a project config file (project.ini)

        :param filepath: Path to write the project config.
        :type filepath: str
        :return: None
        """
        config = configparser.ConfigParser()
        for section, options in self.profile.items():
            config.add_section(section)
            for key, value in options.items():
                config.set(section, key, str(value))
        # Write the configuration to the file
        with open(self.profile["Project"]["project.ini"], 'w') as config_file:
            config.write(config_file)
            print("[Logs]", file=config_file)
            for entry in self.logs:
                print(entry, file=config_file)

    def update_log(self):
        """
        Update the log by adding the new command with a time stamp.

        :return: None
        """
        ctx = click.get_current_context()
        full_command = " ".join(ctx.command_path.split())
        current_datetime = datetime.now()
        formatted_timestamp = current_datetime.strftime('%y%m%d %H:%M')
        new_entry = formatted_timestamp + ": " + full_command
        self.logs.append(new_entry)

    def copy_file(self, source, target):
        """
        Copy a file from source to target with modifying certain key words

        :param source: Path of the source file.
        :type source: str
        :param target: Path of the target file.
        :type target: str
        :return: None
        """
        config_dict = get_dict_from_configs()
        with open(source, 'r') as input_file, open(target, 'w') as output_file:
            for line in input_file:
                # project.ini
                for section, options in self.profile.items():
                    if section != "Logs":
                        for tag, value in options.items():
                            if tag.upper() in line:
                                line = line.replace(tag.upper(), value)
                # configs
                line = replace_variables_by_dict(line, config_dict)
                output_file.write(line)

    def demultiplex(self, method, raw, output):
        """
        Copy the files for demultiplexing according to the given method.

        :param method: One of the methods defined in gpm.config.
        :type method: str
        :param raw: Path of the raw BCL folder.
        :type raw: str
        :param output: Path of the output folder.
        :type output: str
        :return: None
        """
        raw = remove_end_slash(raw)
        raw = os.path.abspath(raw)
        output = remove_end_slash(output)
        output = os.path.abspath(output)
        # Check path for raw data
        if os.path.exists(raw):
            self.profile["Raw data"]["bcl_path"] = raw
        else:
            click.echo("The given path for raw data doesn't exist.")
            click.echo(raw)
            sys.exit()
        # Check path for output path
        raw_name = os.path.basename(raw)
        # If run name is repeated
        if os.path.basename(output) == raw_name:
            click.echo("Please don't repeat the basename of the folder.")
            click.echo("Instead of this:")
            click.echo(output)
            click.echo("Please use this:")
            click.echo(os.path.dirname(output))
            sys.exit()
        # If the run exists, otherwise create a new folder
        if os.path.exists(os.path.join(output, raw_name)):
            click.echo("This run exists in the output directory.")
            click.echo(os.path.join(output, raw_name))
            sys.exit()
        else:
            os.mkdir(os.path.join(output, raw_name))
        # Copy the method
        source_dir = os.path.join(get_gpmdata_path(), "demultiplex", method)
        for filename in os.listdir(source_dir):
            file_path = os.path.join(source_dir, filename)
            target_file = os.path.join(output, raw_name, filename)
            if os.path.isfile(file_path):
                self.copy_file(source=file_path, target=target_file)
        # Update profile
        self.profile["Demultiplexing"]["fastq_path"] = os.path.join(output,
                                                                    raw_name)
        config_path = os.path.join(output, raw_name, PROJECT_INI_FILE)
        self.profile["Project"]["project.ini"] = config_path

    def init_project(self, name):
        """
        Initiate a project by creating a new folder, load the project.ini,
        and populate the files for further processing.

        :param name: Name for the project with the pattern
        YYMMDD_Name1_Name2_Institute_Application
        :type name: str
        :return: None
        """
        # Check name
        check_project_name(name)
        names = name.split("_")
        self.profile["Project"]["date"] = names[0]
        self.profile["Project"]["name1"] = names[1]
        self.profile["Project"]["name2"] = names[2]
        self.profile["Project"]["institute"] = names[3]
        self.profile["Project"]["application"] = names[4]
        # Create project folder
        current_dir = os.getcwd()
        project_path = os.path.join(current_dir, name)
        if os.path.exists(project_path):
            print("The given project path exists already:")
            print(project_path)
        else:
            os.mkdir(project_path)
        # Update project.ini
        self.profile["Project"]["project_path"] = project_path
        config_path = os.path.join(project_path, PROJECT_INI_FILE)
        self.profile["Project"]["project.ini"] = config_path

    def processing(self, method, fastq):
        """
        Copy the files for processing according to the given method.

        :param method: One of the methods defined in gpm.config.
        :type method: str
        :param fastq: Path of the FASTQ folder.
        :type fastq: str
        :return: None
        """
        # Create processing folder
        processing_path = os.path.join(self.profile["Project"]["project_path"],
                                       method)
        if os.path.exists(processing_path):
            click.echo("The folder exists already:")
            click.echo(processing_path)
        else:
            os.mkdir(processing_path)
        # Copy the method
        source_dir = os.path.join(get_gpmdata_path(), "processing", method)
        for filename in os.listdir(source_dir):
            file_path = os.path.join(source_dir, filename)
            target_file = os.path.join(processing_path, filename)
            if os.path.isfile(file_path):
                self.copy_file(source=file_path, target=target_file)
        # Update project.ini
        self.profile["Processing"]["processing_path"] = processing_path
        self.profile["Processing"]["processing_method"] = method
