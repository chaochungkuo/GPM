import os
from os import path
import sys
import click
from collections import OrderedDict
import configparser
from datetime import datetime
from gpm.helper import remove_end_slash, get_gpmdata_path, \
                       check_project_name, get_dict_from_configs, \
                       replace_variables_by_dict, check_analysis_name, \
                       copy_samplesheet
from gpm import PROJECT_INI_FILE

tags_GPM = OrderedDict([("Project", ["date", "name1", "name2", "institute",
                                     "application", "project.ini",
                                     "project_path", "project_name"]),
                        ("Raw data", ["bcl_path"]),
                        ("Demultiplexing", ["demultiplex_path",
                                            "fastq_path",
                                            "fastq_multiqc_path",
                                            "demultiplex_method"]),
                        ("Processing", ["processing_path",
                                        "processing_method",
                                        "processing_qc_path",
                                        "organism",
                                        "genome_assembly"]),
                        ("Analysis", ["analysis_path",
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
            if section == "History":
                logs = config[section]["logs"]
                logs = [log.strip() for log in logs.strip().split("\n")]
                self.logs = logs
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
        self.profile["History"]["logs"] = "\t" + "\n".join(self.logs)
        config = configparser.ConfigParser()
        for section, options in self.profile.items():
            config.add_section(section)
            for key, value in options.items():
                config.set(section, key, str(value))
        # Write the configuration to the file
        with open(self.profile["Project"]["project.ini"], 'w') as config_file:
            config.write(config_file)

    def update_log(self):
        """
        Update the log by adding the new command with a time stamp.

        :return: None
        """
        full_command = sys.argv
        if full_command[0].endswith("/gpm"):
            full_command[0] = "gpm"
        full_command = " ".join(full_command)
        current_datetime = datetime.now()
        formatted_timestamp = current_datetime.strftime('%y%m%d %H:%M')
        new_entry = formatted_timestamp + " >>> " + full_command
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

        :param method: One of the methods defined in gpm.ini.
        :type method: str
        :param raw: Path of the raw BCL folder.
        :type raw: str
        :param output: Path of the output folder.
        :type output: str
        :return: None
        """
        raw = remove_end_slash(raw)
        raw = path.abspath(raw)
        output = remove_end_slash(output)
        output = path.abspath(output)
        # Check path for raw data
        if path.exists(raw):
            self.profile["Raw data"]["bcl_path"] = raw
        else:
            click.echo("The given path for raw data doesn't exist.")
            click.echo(raw)
            sys.exit()
        # Check path for output path
        raw_name = path.basename(raw)
        # If run name is repeated
        if path.basename(output) == raw_name:
            click.echo("Please don't repeat the basename of the folder.")
            click.echo("Instead of this:")
            click.echo(output)
            click.echo("Please use this:")
            click.echo(path.dirname(output))
            sys.exit()
        # If the run exists, otherwise create a new folder
        if path.exists(path.join(output, raw_name)):
            click.echo("This run exists in the output directory.")
            click.echo(path.join(output, raw_name))
            sys.exit()
        else:
            os.mkdir(path.join(output, raw_name))
        # Copy the method
        source_dir = path.join(get_gpmdata_path(), "demultiplex", method)
        for filename in os.listdir(source_dir):
            file_path = path.join(source_dir, filename)
            target_file = path.join(output, raw_name, filename)
            if path.isfile(file_path):
                self.copy_file(source=file_path, target=target_file)
        # Update profile
        demultiplex_path = path.join(output, raw_name)
        self.profile["Demultiplexing"]["demultiplex_path"] = demultiplex_path
        multiqc_path = path.join(demultiplex_path,
                                 "multiqc/multiqc_report.html")
        self.profile["Demultiplexing"]["fastq_multiqc_path"] = multiqc_path
        self.profile["Demultiplexing"]["demultiplex_method"] = method
        config_path = path.join(output, raw_name, PROJECT_INI_FILE)
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
        self.profile["Project"]["project_name"] = name
        self.profile["Project"]["project_string"] = " ".join(names)
        # Create project folder
        current_dir = os.getcwd()
        project_path = path.join(current_dir, name)
        if path.exists(project_path):
            print("The given project path exists already:")
            print(project_path)
        else:
            os.mkdir(project_path)
        # Update project.ini
        self.profile["Project"]["project_path"] = project_path
        config_path = path.join(project_path, PROJECT_INI_FILE)
        self.profile["Project"]["project.ini"] = config_path

    def processing(self, method, fastq):
        """
        Copy the files for processing according to the given method.

        :param method: One of the methods defined in gpm.ini.
        :type method: str
        :param fastq: Path of the FASTQ folder.
        :type fastq: str
        :return: None
        """
        # Create processing folder
        processing_path = path.join(self.profile["Project"]["project_path"],
                                    method)
        if path.exists(processing_path):
            click.echo("The folder exists already:")
            click.echo(processing_path)
        else:
            os.mkdir(processing_path)
        # Copy the method
        source_dir = path.join(get_gpmdata_path(), "processing", method)
        for filename in os.listdir(source_dir):
            file_path = path.join(source_dir, filename)
            target_file = path.join(processing_path, filename)
            if path.isfile(file_path):
                self.copy_file(source=file_path, target=target_file)
        # Update project.ini
        self.profile["Processing"]["processing_path"] = processing_path
        self.profile["Processing"]["processing_method"] = method

    def add_analysis_dir(self):
        """
        Add the analysis folder and update project profile.

        :return: None
        """
        analysis_dir = path.join(self.profile["Project"]["project_path"],
                                 "analysis")
        if not path.exists(analysis_dir):
            os.mkdir(analysis_dir)
        self.profile["Analysis"]["analysis_path"] = analysis_dir

    def add_analysis_report(self, application):
        """
        Add the Rmd report according to application type.

        :return: None
        """
        source_file = path.join(get_gpmdata_path(), "analysis",
                                "Analysis_Report_"+application+".Rmd")
        target_file = path.join(self.profile["Analysis"]["analysis_path"],
                                "Analysis_Report_"+application+".Rmd")
        self.copy_file(source_file, target_file)

    def show_analysis_list(self):
        """
        Show the list of analysis templates available.

        :return: None
        """
        analysis_dict = self.load_analysis_config()
        for group in analysis_dict.keys():
            click.echo(click.style(group, fg='bright_green'))
            for label in analysis_dict[group].keys():
                # click.echo("  <<< " + label + " >>>")
                for i, file in enumerate(analysis_dict[group][label]):
                    if i == 0:
                        click.echo("{:<25} {:<}".format(label,
                                                        file.split("/")[2]))
                    else:
                        click.echo("{:<25} {:<}".format("",
                                                        file.split("/")[2]))
            click.echo("")

    def load_analysis_config(self):
        """
        Load the analysis config file from GPM into a dictionary.

        :return: A dictionary for analysis:file
        """
        analysis_config = path.join(get_gpmdata_path(), "config",
                                    "analysis.config")
        analysis_dict = OrderedDict()
        with open(analysis_config) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                elif len(line.split(",")) == 3:
                    x_list = [x.strip() for x in line.split(",")]
                    if x_list[0] not in analysis_dict:
                        analysis_dict[x_list[0]] = OrderedDict()
                    if x_list[1] not in analysis_dict[x_list[0]]:
                        analysis_dict[x_list[0]][x_list[1]] = []
                    analysis_dict[x_list[0]][x_list[1]].append(x_list[2])
        return analysis_dict

    def add_analysis_template(self, analysis_name):
        """
        Add the files of the defined analysis into analysis folder.

        :return: None
        """
        analysis_dict = self.load_analysis_config()
        check_analysis_name(analysis_dict, analysis_name)
        source_dir = get_gpmdata_path()
        target_dir = self.profile["Analysis"]["analysis_path"]
        for group in analysis_dict.keys():
            for label in analysis_dict[group].keys():
                if label == analysis_name:
                    group_dir = path.join(target_dir, group)
                    if not path.exists(group_dir):
                        os.makedirs(group_dir)
                    for template in analysis_dict[group][label]:
                        click.echo("  "+template)
                        source_file = path.join(source_dir, template)
                        target_file = path.join(group_dir,
                                                path.basename(template))
                        self.copy_file(source_file, target_file)

    def run_analysis_codes(self, analysis_name):
        """
        Run the corresponding codes while inserting new analysis method.

        :return: None
        """
        if analysis_name == "DGEA_RNAseq":
            sheet = path.join(self.profile["Processing"]["processing_path"],
                              "samplesheet.csv")
            dir_analysis = self.profile["Analysis"]["analysis_path"]
            copy_samplesheet(sheet,
                             path.join(dir_analysis, "DGEA",
                                       "samplesheet.csv"))
