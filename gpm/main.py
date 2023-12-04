import click
from gpm.helper import get_gpm_config
from gpm.version import version
from gpm.gpm import GPM

help_messages = {"demultiplex_raw": "Define the folder of BCL files as the "
                 "raw input for demultiplexing.",
                 "demultiplex_output": "Define the output directory where a "
                 "new folder with the same name as raw data will be created."}


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


if __name__ == '__main__':
    main()
