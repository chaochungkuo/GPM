import os
from os import path
import sys
import click
import glob
import shutil
from collections import OrderedDict
import configparser
from datetime import datetime
from gpm.helper import remove_end_slash, get_gpmdata_path, \
                       check_project_name, get_dict_from_configs, \
                       replace_variables_by_dict, check_analysis_name, \
                       copy_samplesheet, get_gpm_config, append_file_to_another
from gpm.messages import show_tree, show_instructions
from gpm import PROJECT_INI_FILE
from gpm.exports import check_export_directory, get_htaccess_path, \
                        htpasswd_create_user
from gpm.project_ini_struc import tags_GPM


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
        filepath = os.path.abspath(filepath)
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
        # Check the project.ini path with symbolic link
        self.prefix = self.symbolic_profile_path(filepath)
        # print("Detected symbolic link:", self.prefix)
        self.exports = self.profile
        if self.prefix != "":
            self.update_with_symlink(self.prefix)

    def symbolic_profile_path(self, filepath):
        if self.profile["Project"]["project.ini"] == filepath:
            return ""
        else:
            if len(filepath) > len(self.profile["Project"]["project.ini"]):
                prefix = filepath.replace(
                    self.profile["Project"]["project.ini"], "")
            else:
                prefix = self.profile["Project"]["project.ini"].replace(
                    filepath, "")
            return prefix

    def update_with_symlink(self, prefix):
        """
        Update the project.ini path with symbolic link.
        """
        for section in self.exports.keys():
            for tag in self.exports[section].keys():
                if self.exports[section][tag].startswith("/"):
                    # self.exports[section][tag] = os.path.join(
                    #     prefix, self.exports[section][tag])
                    self.exports[section][tag] = \
                        prefix + self.exports[section][tag]
                    # print(self.exports[section][tag])

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

    def replace_variables_by_project_ini(self, line):
        for section, options in self.profile.items():
            if section != "Logs":
                for tag, value in options.items():
                    if "PROJECT_"+tag.upper() in line:
                        line = line.replace("PROJECT_"+tag.upper(), value)
        return line

    def replace_variable(self, line, config_dict):
        # project.ini
        line = self.replace_variables_by_project_ini(line)
        # configs
        line = replace_variables_by_dict(line, config_dict)
        return line

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
                line = self.replace_variable(line, config_dict)
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
        # Show instructions
        show_tree(demultiplex_path)
        show_instructions("demultiplex", method)

    def update_project_name(self, name):
        check_project_name(name)
        names = name.split("_")
        self.profile["Project"]["date"] = names[0]
        self.profile["Project"]["name1"] = names[1]
        self.profile["Project"]["name2"] = names[2]
        self.profile["Project"]["institute"] = names[3]
        self.profile["Project"]["application"] = names[4]
        self.profile["Project"]["project_name"] = name
        self.profile["Project"]["project_string"] = " ".join(names)

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
        self.update_project_name(name)
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
        if fastq:
            self.update_fastq_path(fastq)
        # Copy the method
        source_dir = path.join(get_gpmdata_path(), "processing", method)
        for filename in os.listdir(source_dir):
            file_path = path.join(source_dir, filename)
            target_file = path.join(processing_path, filename)
            if path.isfile(file_path):
                self.copy_file(source=file_path, target=target_file)
        # insert nextflow.config
        if "nfcore" in method:
            gpm_nextflow_config = path.join(get_gpmdata_path(),
                                            "config", "nextflow.config")
            process_nextflow_config = path.join(processing_path,
                                                "nextflow.config")
            if path.exists(process_nextflow_config):
                append_file_to_another(gpm_nextflow_config,
                                       process_nextflow_config)
            else:
                self.copy_file(source=gpm_nextflow_config,
                               target=process_nextflow_config)
        # Update project.ini
        self.profile["Processing"]["processing_path"] = processing_path
        self.profile["Processing"]["processing_method"] = method
        # Show instructions
        show_tree(self.profile["Project"]["project_path"])
        show_instructions("processing", method)

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
        files = ["Analysis_Report_"+application+".Rmd",
                 "report_functions.R",
                 "references.bib"]
        for copy_file in files:
            source_file = path.join(get_gpmdata_path(), "analysis",
                                    copy_file)
            target_file = path.join(self.profile["Analysis"]["analysis_path"],
                                    copy_file)
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
                        if os.path.isfile(source_file):
                            self.copy_file(source_file, target_file)
                        elif os.path.isdir(source_file):
                            shutil.copytree(source_dir, target_file)
        # Show instructions
        show_tree(group_dir)
        show_instructions("analysis", analysis_name)

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

    def load_export_config(self):
        self.export_structure = []
        config_dict = get_dict_from_configs()
        cfg_path = os.path.join(get_gpmdata_path(), "config/export.config")
        with open(cfg_path) as config:
            for line in config:
                if line.startswith("#"):
                    continue
                else:
                    ll = [le.strip() for le in line.split(";")]
                    if len(ll) == 4:
                        if (
                            ll[0] == "all" or
                            ll[0].lower() ==
                            self.profile["Project"]["application"].lower()
                        ):
                            ll[1] = self.replace_variable(ll[1], config_dict)
                            self.export_structure.append(ll)

    def export(self, export_dir, tar=False):
        def handle_rename(export_dir, entry):
            # print(os.path.basename(entry[1]))
            if entry[3]:
                target = os.path.join(export_dir, entry[2], entry[3])
            else:
                target = os.path.join(export_dir, entry[2],
                                      os.path.basename(entry[1]))
            return target

        export_dir = os.path.abspath(export_dir)
        check_export_directory(export_dir)
        # Load export name to project.ini
        export_name = os.path.basename(export_dir)
        self.update_project_name(export_name)
        # Creating soft links of the files
        self.load_export_config()
        for entry in self.export_structure:
            # print(entry)
            if not entry[1]:  # make the folder
                target = os.path.join(export_dir, entry[2])
                if not os.path.exists(target):
                    os.makedirs(target)
            else:
                origin_f = os.path.join(
                    self.profile["Project"]["project_path"],
                    entry[1]
                )
                # A directory
                if os.path.isdir(origin_f):
                    target = handle_rename(export_dir, entry)
                    os.symlink(self.prefix+origin_f, target,
                               target_is_directory=True)
                # A file
                elif os.path.isfile(origin_f):
                    target = handle_rename(export_dir, entry)
                    os.symlink(self.prefix+origin_f, target,
                               target_is_directory=False)
                # A pattern for many files
                else:
                    target_dir = os.path.join(export_dir, entry[2])
                    if not os.path.exists(target_dir):
                        os.makedirs(target_dir)
                    for matching_file in glob.glob(origin_f):
                        target = os.path.join(target_dir,
                                              os.path.basename(matching_file))
                        os.symlink(matching_file, target,
                                   target_is_directory=False)

    def add_htaccess(self, export_dir):
        htaccess_path = get_htaccess_path()
        self.copy_file(htaccess_path,
                       os.path.join(export_dir, ".htaccess"))
        # shutil.chown(os.path.join(export_dir, ".htaccess"), group=GROUPNAME)

    def create_user(self, export_dir, raw_export=False):
        export_URL = os.path.join(get_gpm_config("EXPORT", "EXPORT_URL"),
                                  self.profile["Project"]["project_name"])
        export_user, export_password = htpasswd_create_user(
            export_dir,
            export_URL,
            self.profile["Project"]["name1"].lower(),
            self.profile["Project"]["application"])
        self.profile["Export"]["export_URL"] = export_URL
        self.profile["Export"]["export_user"] = export_user
        self.profile["Export"]["export_password"] = export_password

    def update_username(self, username):
        self.profile["Export"]["export_user"] = username

    def update_fastq_path(self, fastq_path):
        self.profile["Demultiplexing"]["fastq_path"] = fastq_path
