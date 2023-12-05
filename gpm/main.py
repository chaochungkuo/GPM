import click
from gpm.helper import get_gpm_config
from gpm.version import version
from gpm.gpm import GPM

help_messages = {"demultiplex_raw": "Define the folder of BCL files as the "
                 "raw input for demultiplexing.",
                 "demultiplex_output": "Define the output directory where a "
                 "new folder with the same name as raw data will be created.",
                 "init_fastq": "Define the path to the FASTQ files for "
                 "initiating a new project.",
                 "init_name": "Define the name of the project in GPM in the "
                 "format of YYMMDD_Name1_Name2_Institute_Application."}


@click.group()
@click.version_option(version)
def main():
    """The Genomic Project Manager is a powerful tool designed
       to streamline and automate bioinformatic workflows within a
       core facility.\n
       Contact: chao-chung.kuo@rwth-aachen.de\n
       Github: https://github.com/chaochungkuo/GPM
    """
    pass


@main.command()
@click.option('-r', '--raw',
              help=help_messages["demultiplex_raw"], required=True)
@click.option('-o', '--output',
              help=help_messages["demultiplex_output"], required=True)
@click.option('-m', '--method',
              type=click.Choice(get_gpm_config("GPM", "DEMULTIPLEX_METHODS"),
                                case_sensitive=False))
def demultiplex(method, raw, output):
    """
    GPM offers various ways for demultiplexing according to different kits or
    sequencers.
    """
    pm = GPM()
    pm.demultiplex(method, raw, output)
    pm.update_log()
    pm.write_project_config_file()


@main.command()
@click.argument('project_config')
@click.option('-fq', '--fastq',
              help=help_messages["init_fastq"], required=True)
@click.option('-n', '--name',
              help=help_messages["init_name"], required=True)
@click.option('-p', '--processing', required=False,
              type=click.Choice(get_gpm_config("GPM", "PROCESSING_METHODS"),
                                case_sensitive=False))
def init(project_config, fastq, name, processing):
    """
    Initiate a project with GPM by inheriting project.ini for further
    processing and analysis.
    """
    pm = GPM()
    pm.load_project_config_file(project_config)
    pm.init_project(name)
    if processing:
        pm.processing(processing, fastq)
    pm.update_log()
    pm.write_project_config_file()


@main.command()
@click.argument('project_config')
@click.option('-fq', '--fastq',
              help=help_messages["init_fastq"], required=True)
@click.option('-p', '--processing',
              type=click.Choice(get_gpm_config("GPM", "PROCESSING_METHODS"),
                                case_sensitive=False))
def processing(project_config, fastq, processing):
    """
    Add additional processing methods in an existed project.
    """
    pm = GPM()
    pm.load_project_config_file(project_config)
    pm.processing(processing, fastq)
    pm.update_log()
    pm.write_project_config_file()


if __name__ == '__main__':
    main()
