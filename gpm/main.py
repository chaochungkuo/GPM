import os
import click
from gpm.archive import archive_folders
from gpm.clean import clean_folders
from gpm.helper import get_gpm_config
from gpm.samplesheet import generate_samples, fastq_dir_to_samplesheet
from gpm.__version__ import __version__
from gpm.gpm import GPM
from gpm.exports import check_export_directory, export_empty_folder, tar_exports

help_messages = {
    "demultiplex_raw": "Define the folder of BCL files as the "
    "raw input for demultiplexing.",
    "demultiplex_output": "Define the output directory where a "
    "new folder with the same name as raw data will be created.",
    "init_fastq": "Define the path to the FASTQ files for " "initiating a new project.",
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
    "(1-based) will be joined to create final sample name.",
    "from-config": "Define the project.ini to inherit from.",
    "authors": "Define the authors of the new project according to gpm.ini."
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
    # os.umask(0o002)
    pass


@main.command()
@click.option("-r", "--raw", help=help_messages["demultiplex_raw"], required=True)
@click.option("-o", "--output", help=help_messages["demultiplex_output"], required=True)
@click.option(
    "-m",
    "--method",
    required=True,
    type=click.Choice(
        get_gpm_config("GPM", "GPM_DEMULTIPLEX_METHODS"), case_sensitive=False
    ),
)
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
@click.option("-fc", "--from-config", required=False, help=help_messages["from-config"])
@click.option("-fq", "--fastq", help=help_messages["init_fastq"], required=False)
@click.option("-n", "--name", help=help_messages["init_name"], required=True)
@click.option(
    "-p",
    "--processing",
    required=False,
    type=click.Choice(
        get_gpm_config("GPM", "GPM_PROCESSING_METHODS"), case_sensitive=False
    ),
    help="Define the pipeline for this project.",
)
@click.option("-a", "--authors", help=help_messages["authors"],
              required=False, default=None)
def init(from_config, fastq, name, processing, authors):
    """
    Initiate a project with GPM by inheriting project.ini for further
    processing and analysis.
    """
    pm = GPM()
    if from_config:
        pm.load_project_config_file(from_config)
    pm.init_project(name, authors)
    if processing:
        pm.processing(processing, fastq)
    if fastq:
        pm.update_fastq_path(fastq)
    pm.update_log()
    pm.write_project_config_file()


@main.command()
@click.argument("project_config")
@click.option("-fq", "--fastq", help=help_messages["init_fastq"], required=True)
@click.option(
    "-p",
    "--processing",
    type=click.Choice(
        get_gpm_config("GPM", "GPM_PROCESSING_METHODS"), case_sensitive=False
    ),
    help="Define the pipeline for this project.",
)
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
@click.argument("project_config")
@click.option(
    "-r",
    "--report",
    type=click.Choice(get_gpm_config("GPM", "GPM_REPORTS"), case_sensitive=False),
    help="Define the kind of report for this project.",
)
@click.option(
    "-ls",
    "--list",
    "show_list",
    required=False,
    default=False,
    is_flag=True,
    help="List all the available analysis templates.",
)
@click.option(
    "-a",
    "--add",
    "add_template",
    required=False,
    default="",
    help="Add the defined analysis template into the project.",
)
def analysis(project_config, report, show_list, add_template):
    """
    Initiate analyses and reports in an existed project.
    """
    pm = GPM()
    pm.load_project_config_file(project_config)
    pm.add_analysis_dir()
    if show_list:
        pm.show_analysis_list()
    else:
        if report:
            pm.add_analysis_report(report)
        if add_template:
            pm.add_analysis_template(add_template)
            # pm.run_analysis_codes(add_template)
        pm.update_log()
        pm.write_project_config_file()


@main.command()
@click.argument("export_folder")
@click.option(
    "-c",
    "--config",
    "config",
    required=False,
    default="",
    help="Define the config file of the project.",
)
# @click.option('-s', '--sym', "prefix", required=False, default="",
#               help="Add the sprefixymbolic prefix to all symbolic paths.")
@click.option(
    "-nc",
    "--no-cloud",
    "no_cloud",
    required=False,
    default=False,
    help="Do not create cloud export (Default: False).",
)
@click.option(
    "-u",
    "--user",
    "username",
    required=False,
    default=None,
    help="Define the user name if needed.",
)
@click.option(
    "-t",
    "--tar",
    "tar",
    required=False,
    default=False,
    is_flag=True,
    help="Tar the folders for download.",
)
@click.option(
    "-g",
    "--gzip",
    "gzip",
    default=False,
    show_default=True,
    is_flag=True,
    help="Generate tar.gz instead of tar.",
)
@click.option(
    "-b",
    "--bcl",
    "bcl",
    required=False,
    default="",
    help="Define the path to BCL folder.",
)
def export(export_folder, config, no_cloud, username, bcl, tar, gzip):
    """
    Export the project to the target folder with symbolic links.
    """

    pm = GPM()
    if config:  
        pm.load_project_config_file(config)
        pm.export(export_folder)
    else:
        pm.export(export_folder, symlink=False)
    if username:
        pm.update_username(username)
    pm.add_htaccess(export_folder)
    pm.create_user(export_folder)  # Here the echo is happening
    if not no_cloud:
        pm.create_cloud_export(os.path.basename(export_folder))
    pm.update_export_report(export_folder)
    pm.echo_export_info()
    pm.echo_wget_msg(export_folder)
    pm.echo_json_info()
    pm.update_log()
    if config:
        pm.write_project_config_file()
    if bcl:
        os.symlink(bcl, os.path.join(export_folder, "BCL"),
                   target_is_directory=True)
    if tar:
        tar_exports(
            export_folder=export_folder, gzip=gzip, dry_run=False, same_server=False
        )


@main.command()
@click.argument("export_folder")
@click.option(
    "-d",
    "--dry-run",
    "dry_run",
    default=False,
    show_default=True,
    is_flag=True,
    help="Dry run without actual execution.",
)
@click.option(
    "-g",
    "--gzip",
    "gzip",
    default=False,
    show_default=True,
    is_flag=True,
    help="Generate tar.gz instead of tar.",
)
def tar_export(export_folder, dry_run, gzip):
    """Tar the sub folders under the export directory with symlinks,
    except compressed_tar folder. This command should be executed on
    export server."""
    tar_exports(
        export_folder=export_folder, gzip=gzip, dry_run=dry_run, same_server=True
    )


@main.command()
@click.argument("target_folders", nargs=-1)
@click.option(
    "-d",
    "--dry-run",
    "dry_run",
    default=False,
    show_default=True,
    is_flag=True,
    help="Dry run without actual execution.",
)
@click.option("-b", "--before-date", "before", default="",
              help="Filter the folders by the date (YYMMDD) in its name. Any folders after this date will be skipped.")
# @click.option("-a", "--after-date", "after", default="",
#               help="Filter the folders by the date in its name.")
@click.option(
    "-v",
    "--v",
    "show_each_file",
    default=False,
    show_default=True,
    is_flag=True,
    help="Show the size of each file.",
)
@click.option(
    "-k",
    "--keep",
    "keep_files",
    default=True,
    show_default=True,
    is_flag=True,
    help="Keep the folders which contain an empty .keep file in the directory.",
)
def clean(target_folders, dry_run, show_each_file, keep_files,before):
    """Clean the given folders by deleting the patterns defined in gpm.ini."""
    if len(target_folders) == 1 and not os.path.isdir(target_folders[0]):
        click.echo("No folders is provided.")
    else:
        click.echo("Following files/folders could be cleaned:")
        clean_folders(target_folders, dry=dry_run,
                      show_each_file=show_each_file,
                      keep_files=keep_files,
                      before=before)


@main.command()
@click.argument("source_folders", nargs=-1)
@click.argument("archive_destination")
@click.option(
    "-d",
    "--dry-run",
    "dry_run",
    default=False,
    show_default=True,
    is_flag=True,
    help="Dry run without actual execution.",
)
@click.option(
    "-v",
    "--verbose",
    "verbose",
    default=False,
    show_default=True,
    is_flag=True,
    help="Show details.",
)
def archive(source_folders, archive_destination, dry_run, verbose):
    """
    Archive the given folders to the archive destination and delete the
    original folders. This should be done after cleaning by gpm clean.
    """
    if len(source_folders) == 1 and not os.path.isdir(source_folders[0]):
        click.echo("No folders is provided.")
    else:
        archive_folders(
            source_folders, archive_destination, dry=dry_run, verbose=verbose
        )


@main.command()
@click.argument("samplesheet")
@click.argument("fastq_dir")
@click.option(
    "-st", default="unstranded", show_default=True, help=help_messages["samplesheet_st"]
)
@click.option(
    "-r1",
    default="_R1_001.fastq.gz",
    show_default=True,
    help="File extension for read 1.",
)
@click.option(
    "-r2",
    default="_R2_001.fastq.gz",
    show_default=True,
    help="File extension for read 2.",
)
@click.option(
    "-se", default=False, show_default=True, help=help_messages["samplesheet_se"]
)
@click.option(
    "-sn", default=False, show_default=True, help=help_messages["samplesheet_sn"]
)
@click.option(
    "-sd",
    default="_",
    show_default=True,
    help="Delimiter to use to sanitise sample name.",
)
@click.option("-si", default=1, show_default=True, help=help_messages["samplesheet_si"])
def samplesheet_rnaseq(samplesheet, fastq_dir, st, r1, r2, se, sn, sd, si):
    """Generate sample sheet for nf-core RNAseq pipeline."""
    generate_samples(st, fastq_dir, samplesheet, r1, r2, se, sn, sd, si)


@main.command()
@click.argument("samplesheet")
@click.argument("fastq_dir")
@click.option(
    "-r1",
    default="_R1_001.fastq.gz",
    show_default=True,
    help="File extension for read 1.",
)
@click.option(
    "-r2",
    default="_R2_001.fastq.gz",
    show_default=True,
    help="File extension for read 2.",
)
@click.option(
    "-se", default=False, show_default=True, help=help_messages["samplesheet_se"]
)
@click.option(
    "-sn", default=False, show_default=True, help=help_messages["samplesheet_sn"]
)
@click.option(
    "-sd",
    default="_",
    show_default=True,
    help="Delimiter to use to sanitise sample name.",
)
@click.option("-si", default=1, show_default=True, help=help_messages["samplesheet_si"])
def samplesheet_scrnaseq(samplesheet, fastq_dir, r1, r2, se, sn, sd, si):
    """Generate sample sheet for nf-core single cell RNAseq pipeline."""
    fastq_dir_to_samplesheet(
        fastq_dir=fastq_dir,
        samplesheet_file=samplesheet,
        strandedness="",
        read1_extension=r1,
        read2_extension=r2,
        single_end=se,
        sanitise_name=sn,
        sanitise_name_delimiter=sd,
        sanitise_name_index=si,
        sc=True,
    )


if __name__ == "__main__":
    main()
