import os
import re
import click
from gpm.helper import get_gpm_config


def clean_folders(target_folders, show_total, show_each_file):
    pattern_strings = get_gpm_config("CLEAN", "PATTERNS")
    regex_patterns = [re.compile(pattern) for pattern in pattern_strings]

    print(regex_patterns)
    for folder in target_folders:
        if os.path.isdir(folder):
            matching_files = search_files_by_patterns(folder,
                                                      regex_patterns)
        if matching_files:
            click.echo(click.style(folder, fg='bright_green'))
            if show_each_file:
                for matching_file in matching_files:
                    click.echo(matching_file)
            if show_total:
                total_size = get_total_size(matching_files)
                click.echo(total_size)


def search_files_by_patterns(root_folder, patterns):
    matching_files = []
    for root, dirs, files in os.walk(root_folder):
        matching_files.extend(
            os.path.join(root, file)
            for file in dirs+files
            if any(re.match(re.escape(pattern), file) for pattern in patterns)
        )
    return matching_files


def get_total_size(files):
    total_size = sum(os.path.getsize(file) for file in files)
    return total_size
