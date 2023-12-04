import os
import sys
import click
from collections import OrderedDict
import configparser
from datetime import datetime
from gpm.helper import remove_end_slash
from gpm import PROJECT_INI_FILE

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

    def write_project_config_file(self, filepath):
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
        with open(filepath, 'w') as config_file:
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
        print(source)
        print(target)
        with open(source, 'r') as input_file, open(target, 'w') as output_file:
            for line in input_file:
                for section, options in self.profile.items():
                    if section != "Logs":
                        for tag, value in options.items():
                            if tag.upper() in line:
                                line = line.replace(tag.upper(), value)
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
        source_dir = os.path.dirname(__file__)
        source_dir = os.path.join(source_dir, "data/demultiplex", method)
        for filename in os.listdir(source_dir):
            file_path = os.path.join(source_dir, filename)
            target_file = os.path.join(output, raw_name, filename)
            if os.path.isfile(file_path):
                self.copy_file(source=file_path, target=target_file)
        # Update profile
        self.profile["Demultiplexing"]["fastq_path"] = output
        # Update log
        self.update_log()
        # Write project.ini
        config_path = os.path.join(output, raw_name, PROJECT_INI_FILE)
        self.write_project_config_file(config_path)
