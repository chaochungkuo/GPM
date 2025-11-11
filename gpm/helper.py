import configparser
import datetime
import http
import os
import sys
from os import environ, getenv, path
from urllib.parse import urljoin, urlparse

import click
import pandas as pd
import requests
from bs4 import BeautifulSoup
from requests.auth import HTTPBasicAuth

from gpm import CONFIG_LIST


def get_gpmdata_path():
    """
    Get the GPMDATA path from the environment.
    """
    environ["GPMDATA"] = "/home/mmabrouk/Documents/Projects/gpmdata/"
    if environ.get("GPMDATA"):
        gpm_data_location = path.expanduser(getenv("GPMDATA"))
    else:
        gpm_data_location = path.expanduser(path.join(getenv("HOME"), "gpmdata"))
    return gpm_data_location


def check_project_name(name):
    """
    Check the name of the project with the pattern
    YYMMDD_Name1_Name2_Institute_Application
    """
    split_name = name.split("_")
    if len(split_name) != 5:
        print("Error: Please follow the pattern below:")
        print("YYMMDD_Name1_Name2_Institute_Application")
        sys.exit()
    # date
    sequencing_date = split_name[0]
    try:
        datetime.datetime.strptime(sequencing_date, "%y%m%d")
    except ValueError:
        click.echo("Incorrect date string format. It should be YYMMDD")
        sys.exit()
    # app
    app = split_name[4]
    app_options = [x.lower() for x in get_gpm_config("GPM", "GPM_APPLICATIONS")]
    if app.lower() not in app_options:
        click.echo("Unsupported application. Please take one from below:")
        click.echo(", ".join(get_gpm_config("GPM", "GPM_APPLICATIONS")))
        sys.exit()


def remove_end_slash(path):
    """
    Remove the ending slash in the given path.

    :param path: A path.
    :type path: str
    :return: the modified path
    :rtype: str
    """
    if path.endswith("/"):
        path = path.rstrip("/")
    return path


#################################################################
# Config relevant
#################################################################
def get_config_values(config_name):
    """
    Return a dictionary of the key value pairs in the defined config file.
    User-defined config (*.ini.user) has a higher priority
    than default one (*.ini).
    """
    config_path = path.join(get_gpmdata_path(), "config/" + config_name + ".user")
    if not path.exists(config_path):
        config_path = path.join(get_gpmdata_path(), "config/" + config_name)
    config = configparser.ConfigParser()
    config.read(config_path)
    combined_dict = {}
    for section_name in config.sections():
        section_dict = dict(config.items(section_name))
        combined_dict.update(section_dict)
    return combined_dict


def get_config_section(config_name, section):
    """
    Return a dictionary of the key value pairs in a section of the defined config file.
    User-defined config (*.ini.user) has a higher priority
    than default one (*.ini).
    """
    config_path = path.join(get_gpmdata_path(), "config/" + config_name + ".user")
    if not path.exists(config_path):
        config_path = path.join(get_gpmdata_path(), "config/" + config_name)
    config = configparser.ConfigParser()
    config.read(config_path)
    section_dict = dict(config.items(section))
    return section_dict


def get_config_value(config_name, section, item):
    """
    Return value of the defined item in the defined config file.
    User-defined config (*.ini.user) has a higher priority
    than default one (*.ini).
    """
    config_path = path.join(get_gpmdata_path(), "config/" + config_name + ".user")
    if not path.exists(config_path):
        config_path = path.join(get_gpmdata_path(), "config/" + config_name)
    config = configparser.ConfigParser()
    config.read(config_path)
    res = config[section][item]
    if "," in res:
        res = [x.strip() for x in res.split(",")]
    return res


def get_gpm_config(section, item):
    """
    Return the config from GPMDATA/gpm.ini as a dictionary.
    User-defined config (gpm.ini.user) has a higher priority
    than default one (gpm.ini).
    """
    config_name = "gpm.ini"
    res = get_config_value(config_name, section, item)
    return res


def get_environment_config(section, item):
    """
    Return the config from GPMDATA/environment.ini as a dictionary.
    User-defined config (environment.ini.user) has a higher priority
    than default one (environment.ini).
    """
    config_name = "environment.ini"
    res = get_config_value(config_name, section, item)
    return res


def get_dict_from_configs():
    combined_dict = {}
    for sel_config in CONFIG_LIST:
        config_dict = get_config_values(sel_config)
        combined_dict.update(config_dict)
    # print(combined_dict)
    return combined_dict


def replace_variables_by_dict(line, input_dict):
    for key, value in input_dict.items():
        if key.upper() in line:
            line = line.replace(key.upper(), value)
    return line


def check_analysis_name(analysis_dict, analysis_name):
    all_names = []
    for k, g in analysis_dict.items():
        all_names = all_names + list(g.keys())
    if analysis_name not in all_names:
        click.echo("Please choose an analysis from the list below")
        click.echo(all_names)
        sys.exit()


def copy_samplesheet(source_samplesheet, target_samplesheet):
    try:
        df = pd.read_csv(source_samplesheet)
        df.drop(columns=["fastq1", "fastq2", "strandness"], inplace=True)
        df[["label1", "label2"]] = df["sample"].str.split("_", expand=True)
        df.to_csv(target_samplesheet, index=False)
    except Exception:
        warning = (
            "No samplesheet from processing path is available.\n"
            "Please generate the samplesheet.csv manually."
        )
        click.echo(click.style(warning, fg="red"))


def append_file_to_another(file1, file2):
    with open(file1, "r") as source_file:
        source_content = source_file.read()
    # Append the content to the destination file
    with open(file2, "a") as destination_file:
        destination_file.write(source_content)


def get_authors(short_names):
    gpm_authors = get_config_section("gpm.ini", "AUTHORS")
    authors = list(gpm_authors.keys())
    res = []
    if short_names is not None:  # authors are defined
        list_short_names = short_names.split(",")
        for name in list_short_names:
            if name in authors:
                res.append(gpm_authors[name])
            else:
                print(f"{name} is not defined in gpm.ini. Skipped.")
    else:  # Not defined and take all available authors
        for name in authors:
            res.append(gpm_authors[name])
    return res


def author_list2string(authors_list, format):
    if format == "RMD":
        authors = ""
        for au in authors_list:
            authors += "  - " + au + "\n"
    elif format == "ipynb":
        authors = []
        for au in authors_list:
            authors.append("  - " + au)
        authors = '\\n",\n    "'.join(authors)
    return authors


def check_links_in_html(html_content, base_url=None):
    """
    Check all links (relative or absolute) in the given HTML content.

    Args:
        html_content (str): HTML text as a string.
        base_url (str, optional): The base URL or file path to resolve relative links.

    Returns:
        list of dict: List containing link info with status.
    """
    soup = BeautifulSoup(html_content, "html.parser")
    links = [a.get("href") for a in soup.find_all("a", href=True)]

    results = []

    for link in links:
        if not link:
            continue

        # Handle relative URLs
        full_url = urljoin(base_url, link) if base_url else link

        # Check if it's a local file
        parsed = urlparse(full_url)
        if parsed.scheme in ("http", "https"):
            try:
                response = requests.head(full_url, allow_redirects=True, timeout=5)
                status = (
                    "OK"
                    if response.status_code < 400
                    else f"Error {response.status_code}"
                )
            except Exception as e:
                status = f"Failed ({e})"
        else:
            # Assume it's a local file path
            local_path = parsed.path
            if os.path.exists(local_path):
                status = "OK"
            else:
                status = "File Not Found"

        results.append(
            {
                "link": link,
                "resolved_link": full_url,
                "status": status,
            }
        )

    return results


def find_all_linked_htmls(start_html, visited=None):
    """
    Recursively find all HTML files linked from the start HTML file,
    handling relative paths properly.

    Args:
        start_html (str): Starting HTML file path or URL.
        visited (set, optional): Set of visited HTML files.

    Returns:
        set: All reachable HTML file paths or URLs.
    """
    if visited is None:
        visited = set()

    if start_html in visited:
        return visited

    print(f"Visiting: {start_html}")
    visited.add(start_html)

    parsed = urlparse(start_html)

    # Load HTML content
    try:
        if parsed.scheme in ("http", "https"):
            response = requests.get(start_html, timeout=5)
            response.raise_for_status()
            html_content = response.text
            base_for_links = start_html  # URL base
        else:
            local_path = parsed.path
            with open(local_path, "r", encoding="utf-8") as f:
                html_content = f.read()
            base_for_links = (
                "file://" + os.path.dirname(os.path.abspath(local_path)) + "/"
            )  # Local base
    except Exception as e:
        print(f"Failed to open {start_html}: {e}")
        return visited

    # Parse and find links
    soup = BeautifulSoup(html_content, "html.parser")
    links = [a.get("href") for a in soup.find_all("a", href=True)]

    for link in links:
        if not link:
            continue

        # Resolve link relative to the current page
        resolved_link = urljoin(base_for_links, link)

        # Only process HTML files
        if resolved_link.endswith((".html", ".htm")):
            if resolved_link not in visited:
                find_all_linked_htmls(resolved_link, visited=visited)

    return visited


def get_api_creds() -> HTTPBasicAuth | None:
    api_pass = os.getenv("GPM_PASS")
    return HTTPBasicAuth("GPM", api_pass) if api_pass else None


def query_api(
    endpoint: str,
) -> requests.Response:
    """
    Query the GPM API with the given URL and parameters.
    """
    auth = get_api_creds()
    response = requests.get(endpoint, auth=auth)
    if response.status_code != http.HTTPStatus.OK:
        raise Exception(
            f"API request failed with status code {response.status_code} for {endpoint} and {auth}"
        )
    return response


def get_flowcell_id(dir: str) -> str:
    return os.path.basename(dir).split("_")[-1][1:]
