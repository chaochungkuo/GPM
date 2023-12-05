import click
from gpm.helper import get_gpm_config, generate_samples
from gpm.__version__ import __version__
from gpm.gpm import GPM

help_messages = {"demultiplex_raw": "Define the folder of BCL files as the "
                 "raw input for demultiplexing.",
                 "demultiplex_output": "Define the output directory where a "
                 "new folder with the same name as raw data will be created.",
                 "init_fastq": "Define the path to the FASTQ files for "
                 "initiating a new project.",
                 "init_name": "Define the name of the project in GPM in the "
                 "format of YYMMDD_Name1_Name2_Institute_Application.",
                 "samplesheet_st": "Value for 'strandedness' in samplesheet. "
                 "Must be one of 'unstranded', 'forward', 'reverse'.",
                 "samplesheet_se": "Single-end information will be "
                 "auto-detected but this option forces paired-end FastQ files "
                 "to be treated as single-end so only read 1 information is "
                 "included in the samplesheet.",
                 "samplesheet_sn": "Whether to further sanitise FastQ file "
                 "name to get sample id. Used in conjunction with "
                 "--sanitise_name_delimiter and --sanitise_name_index.",
                 "samplesheet_si": "After splitting FastQ file name by "
                 "--sanitise_name_delimiter all elements before this index "
                 "(1-based) will be joined to create final sample name."
                 }


@click.group()
@click.version_option(__version__)
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


@main.command()
@click.argument('samplesheet')
@click.argument('fastq_dir')
@click.option('-st', default="unstranded", show_default=True,
              help=help_messages["samplesheet_st"])
@click.option('-r1', default="_R1_001.fastq.gz", show_default=True,
              help="File extension for read 1.")
@click.option('-r2', default="_R2_001.fastq.gz", show_default=True,
              help="File extension for read 2.")
@click.option('-se', default=False, show_default=True,
              help=help_messages["samplesheet_se"])
@click.option('-sn', default=False, show_default=True,
              help=help_messages["samplesheet_sn"])
@click.option('-sd', default="_", show_default=True,
              help="Delimiter to use to sanitise sample name.")
@click.option('-si', default=1, show_default=True,
              help=help_messages["samplesheet_si"])
def samplesheet_rnaseq(samplesheet, fastq_dir, st, r1, r2, se, sn, sd, si):
    """Generate sample sheet for nf-core RNAseq pipeline."""
    generate_samples(st, fastq_dir, samplesheet,
                     r1, r2, se, sn, sd, si)


if __name__ == '__main__':
    main()
