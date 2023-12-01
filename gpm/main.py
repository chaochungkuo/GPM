import click
from gpm.version import version
help_messages = {"": ""}

###################################################################
# Main function
###################################################################


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
