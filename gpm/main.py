import click
import os
import sys
from gpm.helper import load_gpm_config
from gpm.version import version

help_messages = {"demultiplex_raw": "Define the folder of BCL files as the "
                 "raw input for demultiplexing.",
                 "demultiplex_output": "Define the output directory where a "
                 "new folder with the same name as raw data will be created."}

gpm_config = load_gpm_config()
print(gpm_config["GPM"]["DEMULTIPLEX_METHODS"])


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
              type=click.Choice(gpm_config["GPM"]["DEMULTIPLEX_METHODS"],
                                case_sensitive=False))
def demultiplex(method, raw, output):
    """
    GPM offers various ways for demultiplexing according to different kits or
    sequencers.
    """
    # Check output directory
    if not os.path.exists(output):
        click.echo("Output directory doesn't exist.")
        sys.exit()
    # Check existed folder
    raw_name = os.path.basename(raw)
    if os.path.exists(os.path.join(output, raw_name)):
        click.echo("This run exists in the output directory.")
        click.echo(os.path.join(output, raw_name))
        sys.exit()

    # TODO: generate
