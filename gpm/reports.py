import os
import subprocess
import click
from bs4 import BeautifulSoup
from typing import Dict, List, Set
from urllib.parse import urlparse, urljoin
import requests

def extract_html_titles(root_path: str) -> Dict[str, List[str]]:
    """
    Recursively processes HTML files under the given path and extracts their titles.
    
    Args:
        root_path (str): The root directory to start searching from
        
    Returns:
        Dict[str, List[str]]: A dictionary where:
            - key: HTML filename
            - value: [relative_path, heading_title]
    """
    result = {}
    
    for root, _, files in os.walk(root_path):
        for file in files:
            if file.endswith('.html'):
                file_path = os.path.join(root, file)
                relative_path = os.path.relpath(file_path, root_path)
                
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        soup = BeautifulSoup(f.read(), 'html.parser')
                        title = soup.title.string if soup.title else None
                        
                        if not title:
                            # Use filename without extension and replace underscores with spaces
                            title = os.path.splitext(file)[0].replace('_', ' ')
                        
                        result[file] = [relative_path, title]
                except Exception as e:
                    print(f"Error processing {file_path}: {str(e)}")
                    result[file] = [relative_path, "Error reading file"]
    
    return result

def inserting_sub_reports(analysis_dir: str, report: str, sub_reports: Dict[str, List[str]]) -> None:
    """
    Process an Rmd file to uncomment lines containing sub-report keys and insert remaining reports.
    
    Args:
        analysis_dir (str): Directory containing the Rmd file
        report (str): Name of the report (without extension)
        sub_reports (Dict[str, List[str]]): Dictionary of sub-reports where:
            - key: HTML filename
            - value: [relative_path, title]
    """
    rmd_file = os.path.join(analysis_dir, f"Analysis_Report_{report}.Rmd")
    temp_file = rmd_file + ".tmp"
    
    try:
        with open(rmd_file, 'r', encoding='utf-8') as f_in, open(temp_file, 'w', encoding='utf-8') as f_out:
            found_reports = []
            for line in f_in:
                # Check if line contains any sub-report key
                for key in sub_reports:
                    if key in line:
                        found_reports.append(key)
                        # Uncomment the line by removing HTML comments
                        line = line.replace('<!--', '').replace('-->', '')
                        line = line.lstrip()
                        break
                
                # Check for insertion point
                if '<!-- inserting points -->' in line:
                    f_out.write(line)
                    # Add remaining sub-reports
                    for key, (relative_path, title) in sub_reports.items():
                        if key not in found_reports:
                            f_out.write(f"### [{title}]({relative_path})\n")
                else:
                    f_out.write(line)
        
        # Replace original file with updated version
        os.replace(temp_file, rmd_file)
        
    except Exception as e:
        print(f"Error processing Rmd file: {str(e)}")
        if os.path.exists(temp_file):
            os.remove(temp_file)

def render_rmd_to_html(rmd_file: str) -> None:
    """
    Render an R Markdown file to HTML using R's rmarkdown package.
    
    Args:
        rmd_file (str): Path to the R Markdown file to render
    """
    try:
        command = (
            "source /opt/miniforge3/etc/profile.d/conda.sh && conda activate R4.3 && "
            f"Rscript -e \"rmarkdown::render('{rmd_file}', output_format = 'html_document')\" && "
            "conda deactivate"
        )
        subprocess.run(command, shell=True, executable="/bin/bash")
    except Exception as e:
        print(f"Error rendering R Markdown file: {str(e)}")

def validate_html_links(start_html: str) -> None:
    """
    Recursively validate all HTML links in a file and its linked local HTML files.
    
    Args:
        start_html (str): Path to the starting HTML file
    """
    processed_files = set()
    invalid_links = set()
    
    def process_file(html_path: str) -> None:
        """Process a single HTML file and its links."""
        if html_path in processed_files:
            return
            
        processed_files.add(html_path)
        
        try:
            with open(html_path, 'r', encoding='utf-8') as f:
                soup = BeautifulSoup(f.read(), 'html.parser')
                
                # Get the directory of the current file for resolving relative paths
                current_dir = os.path.dirname(os.path.abspath(html_path))
                
                # Process all links
                for link in soup.find_all('a'):
                    href = link.get('href')
                    if not href:
                        continue
                        
                    # Skip anchor links and external URLs
                    if href.startswith('#') or href.startswith(('http://', 'https://')):
                        continue
                        
                    # Resolve relative path
                    absolute_path = os.path.abspath(os.path.join(current_dir, href))
                    
                    # Check if it's an HTML file
                    if not href.lower().endswith('.html'):
                        continue
                        
                    # Check if file exists
                    if not os.path.exists(absolute_path):
                        invalid_links.add((html_path, href))
                        continue
                        
                    # Recursively process the linked file
                    process_file(absolute_path)
                    
        except Exception as e:
            click.secho(f"Error processing {html_path}: {str(e)}", fg='red')
    
    # Start processing from the initial file
    process_file(os.path.abspath(start_html))
    
    # Report invalid links
    if invalid_links:
        click.secho("\nInvalid links found:", fg='red', bold=True)
        for file_path, link in invalid_links:
            click.secho(f"File: {file_path}", fg='red')
            click.secho(f"Invalid link: {link}\n", fg='red')
    else:
        click.secho("\nAll links are valid!", fg='green', bold=True)
