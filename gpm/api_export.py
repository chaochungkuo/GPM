import glob
import json
import os
import subprocess
import sys
import time
from typing import Any, Dict, Optional, Tuple

import click
import requests

from gpm.exports import generate_password
from gpm.helper import get_gpm_config


def _determine_host():
    """
    Determine the host name from the system using the `hostname` command in Linux.

    Returns:
        str: Host name (e.g., "nextgen", "nextgen2", "nextgen3")
    """
    try:
        result = subprocess.run(
            ["hostname"], capture_output=True, text=True, check=True
        )
        hostname = result.stdout.strip()
        # Extract hostname if it contains dots
        if "." in hostname:
            hostname = hostname.split(".")[0]
        return hostname
    except (subprocess.CalledProcessError, FileNotFoundError):
        click.echo(
            click.style(
                "Unable to determine the host name using the `hostname` command",
                fg="red",
            )
        )
        sys.exit(1)


def _determine_report_section(dest_path: str) -> str:
    """
    Determine the appropriate report_section based on the destination path.

    Maps destination paths to ReportSection enum values:
    - Paths containing '1_Raw_data' → 'raw'
    - Paths containing '2_Processed_data' → 'processed'
    - Paths containing '3_Reports' → 'reports'
    - All other paths → 'general' (default)

    Args:
        dest_path: Destination path (string)

    Returns:
        str: Report section value ('raw', 'processed', 'reports', or 'general')
    """
    dest_path_lower = dest_path.lower()

    if "1_raw_data" in dest_path_lower:
        return "raw"
    elif "2_processed_data" in dest_path_lower:
        return "processed"
    elif "3_reports" in dest_path_lower:
        return "reports"
    else:
        return "general"


def _extract_username_from_project(project_name):
    """
    Extract username from project name.
    Expected format: date_username_group_institute_application
    Returns:
        str: Extracted username in lowercase
    """
    project_parts = project_name.split("_")
    if len(project_parts) >= 2:
        return project_parts[1].lower()
    return project_name.lower()


def _get_backend_config():
    try:
        backends_str = get_gpm_config("EXPORT_ENGINE", "EXPORT_ENGINE_BACKENDS")
        if isinstance(backends_str, str):
            return [b.strip() for b in backends_str.split(",")]
        elif isinstance(backends_str, list):
            return backends_str
    except Exception:
        pass
    return ["apache", "owncloud"]


def extract_authors_from_profile(profile):
    authors = []
    if "authors" in profile.get("Project", {}) and profile["Project"]["authors"]:
        authors_list = profile["Project"]["authors"]
        if isinstance(authors_list, list):
            authors = [str(a) for a in authors_list]
        elif isinstance(authors_list, str):
            authors = [authors_list]
    return authors


def _process_export_entry(entry, project_path, prefix, host, project_name):
    """
    Process a single export structure entry.

    Args:
        entry: Export entry [app, source_path, target_dir, rename]
        project_path: Base project path
        prefix: Symbolic link prefix
        host: Host identifier
        project_name: Project name

    Returns:
        list: List of FileExport entries (may be multiple for glob patterns)
    """
    # Skip entries without source (folder creation only)
    if not entry[1]:
        return []

    source_path = entry[1]
    target_dir = entry[2]
    rename = entry[3] if len(entry) > 3 else None

    # Build absolute source path
    if os.path.isabs(source_path):
        abs_source = prefix + source_path
    else:
        abs_source = prefix + os.path.join(project_path, source_path)

    # Handle glob patterns
    if "*" in source_path or "?" in source_path:
        raise ValueError("Glob patterns are not supported for API export")

    # Single file or directory
    dest_path = os.path.join(
        target_dir, rename if rename else os.path.basename(abs_source)
    )

    section = _determine_report_section(dest_path) if entry[0] != "FASTQ" else "raw"

    return [
        {
            "src": abs_source,
            "dest": dest_path,
            "host": host,
            "project": project_name,
            "mode": "symlink",
            "include_in_report": True,
            "report_section": section,
        }
    ]


def _build_export_list(export_structure, project_path, prefix, host, project_name):
    """
    Convert export structure to list of FileExport entries.

    Args:
        export_structure: List of export entries
        project_path: Base project path
        prefix: Symbolic link prefix
        host: Host identifier
        project_name: Project name

    Returns:
        list: List of FileExport dictionaries
    """
    export_list = []
    for entry in export_structure:
        entries = _process_export_entry(entry, project_path, prefix, host, project_name)
        export_list.extend(entries)
    return export_list


def convert_export_structure_to_job_spec(export_structure, profile, prefix=""):
    """
    Convert GPM's export structure to export_engine's ExportJobSpec format.

    Args:
        export_structure: List of export entries, each with [app, source, target_dir, rename]
        profile: Project profile dictionary
        prefix: Symbolic link prefix (for handling symlinks)

    Returns:
        dict: ExportJobSpec in the format expected by the export engine API
    """
    # Extract project information
    project_name = profile["Project"]["project_name"]
    project_path = profile["Project"]["project_path"]

    # Generate credentials and get configuration
    username = _extract_username_from_project(project_name)
    password = generate_password()
    backend = _get_backend_config()
    authors = extract_authors_from_profile(profile)
    host = _determine_host()

    # Build export list
    export_list = _build_export_list(
        export_structure, project_path, prefix, host, project_name
    )

    # Create ExportJobSpec dict
    job_spec = {
        "project_name": project_name,
        "export_list": export_list,
        "backend": backend,
        "username": username,
        "password": password,
        "authors": authors if authors else None,
    }

    return job_spec


def submit_export_to_api(job_spec: Dict[str, Any]) -> Tuple[Optional[str], Any]:
    """
    Submit export job specification to the export engine API.

    Args:
        job_spec: ExportJobSpec dictionary

    Returns:
        tuple: (job_id, status_dict) or (None, error_message) on failure
    """
    try:
        api_url = get_gpm_config("EXPORT_ENGINE", "EXPORT_ENGINE_API_URL")
        if not api_url:
            return None, "EXPORT_ENGINE_API_URL not configured in gpm.ini"

        # Ensure URL doesn't end with /
        if isinstance(api_url, str):
            api_url = api_url.rstrip("/")
        else:
            return None, "EXPORT_ENGINE_API_URL must be a string"
        endpoint = f"{api_url}/export"

        # Make POST request with timeout
        try:
            response = requests.post(endpoint, json=job_spec, timeout=30)
        except requests.exceptions.Timeout:
            return None, "Request to export engine API timed out after 30 seconds"
        except requests.exceptions.ConnectionError as e:
            return None, f"Failed to connect to export engine API: {str(e)}"
        except requests.exceptions.RequestException as e:
            return None, f"Request to export engine API failed: {str(e)}"

        # Check response
        if response.status_code in (200, 201):
            try:
                result = response.json()
                job_id = result.get("job_id")
                if not job_id:
                    return None, "API response missing job_id"
                status = result.get("status", "submitted")
                return job_id, {"status": status, "response": result}
            except (ValueError, json.JSONDecodeError):
                # Response is not JSON
                return None, f"API returned non-JSON response: {response.text[:200]}"
        else:
            error_msg = f"API request failed with status {response.status_code}"
            try:
                error_detail = response.json()
                error_msg += f": {error_detail}"
            except (ValueError, json.JSONDecodeError):
                error_msg += f": {response.text[:200]}"
            return None, error_msg

    except Exception as e:
        return None, f"Unexpected error submitting to API: {str(e)}"


def poll_final_message(job_id: str, api_url: str) -> Optional[Dict[str, Any]]:
    """
    Poll the final message endpoint for job completion notification.

    Polls the /export/final_message/{job_id} endpoint every 2 seconds
    for up to 60 seconds to retrieve the completion notification.

    Args:
        job_id: The job ID to poll for
        api_url: Base API URL

    Returns:
        Completion notification dictionary if available, None if timeout or error
    """
    api_url = api_url.rstrip("/")
    endpoint = f"{api_url}/export/final_message/{job_id}"

    # Poll for maximum 60 seconds, every 2 seconds
    max_duration = 120
    poll_interval = 2
    start_time = time.time()

    click.echo(
        click.style(
            f"Polling for final message (up to {max_duration} seconds)...",
            fg="bright_blue",
        )
    )

    while time.time() - start_time < max_duration:
        try:
            response = requests.get(endpoint, timeout=5)

            if response.status_code == 200:
                # Success - return the notification
                try:
                    notification = response.json()
                    click.echo(
                        click.style(
                            "Final message retrieved successfully!", fg="bright_green"
                        )
                    )
                    status = notification.get("status", "")
                    msg = notification.get("message", "")
                    formatted_msg = notification.get("formatted_message", "")

                    # Display messages to user
                    if status == "completed":
                        # For completion messages, display both formatted_message and message
                        if formatted_msg:
                            click.echo(click.style(formatted_msg, fg="bright_blue"))
                        if msg:
                            click.echo(
                                click.style(
                                    "\nPlease use the following information for submitting in MS Planner:",
                                    fg="bright_green",
                                )
                            )
                            click.echo(
                                click.style(
                                    msg,
                                    fg="bright_green",
                                )
                            )

                    return notification
                except (ValueError, json.JSONDecodeError):
                    click.echo(
                        click.style(
                            "Warning: Received non-JSON response from final message endpoint",
                            fg="yellow",
                        )
                    )
                    return None

            elif response.status_code == 404:
                # Job ID does not exist
                click.echo(
                    click.style(f"Error: Job ID '{job_id}' does not exist", fg="red")
                )
                return None

            elif response.status_code == 425:
                # Job has not completed yet - continue polling
                time.sleep(poll_interval)
                continue

            else:
                # Other error
                click.echo(
                    click.style(
                        f"Error polling final message: HTTP {response.status_code}",
                        fg="red",
                    )
                )
                return None

        except requests.exceptions.Timeout:
            # Request timeout - continue polling
            time.sleep(poll_interval)
            continue
        except requests.exceptions.ConnectionError as e:
            click.echo(
                click.style(f"Connection error while polling: {str(e)}", fg="red")
            )
            return None
        except requests.exceptions.RequestException as e:
            click.echo(click.style(f"Request error while polling: {str(e)}", fg="red"))
            return None

    # Timeout reached
    click.echo(
        click.style(f"Polling timeout after {max_duration} seconds", fg="yellow")
    )
    return None


def monitor_job_via_polling(
    job_id: str, api_url: str, poll_interval: float = 2.0
) -> dict[str, Any] | None:
    pass


def extract_credentials_from_completion(
    completion_data: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Extract credentials and URLs from completion notification.

    The completion notification contains structured fields: username, password,
    main_report, and publishingResults.

    Args:
        completion_data: Dictionary containing completion notification

    Returns:
        Dictionary with extracted credentials: username, password, export_URL, report_URL, download_url
    """
    if not completion_data or "notification" not in completion_data:
        return {}

    notification = completion_data["notification"]
    credentials = {}

    # Extract from structured fields (always present)
    if "username" in notification:
        credentials["username"] = notification["username"]
    if "password" in notification:
        credentials["password"] = notification["password"]
    if "main_report" in notification:
        credentials["report_URL"] = notification["main_report"]

    # Extract additional data from publishingResults
    publishing_results = notification.get("publishingResults", [])
    if isinstance(publishing_results, list):
        for result in publishing_results:
            # Extract export URLs and download URLs from publishing results
            if isinstance(result, dict):
                # Look for URL fields in publishing results
                if "url" in result:
                    # Determine if it's an export URL or download URL based on backend
                    backend = result.get("backend", "").lower()
                    if backend == "apache":
                        credentials["export_URL"] = result["url"]
                    elif backend in ["owncloud", "cloud"]:
                        credentials["download_url"] = result["url"]
                # Also check for other common field names
                if "export_url" in result:
                    credentials["export_URL"] = result["export_url"]
                if "download_url" in result:
                    credentials["download_url"] = result["download_url"]
                if "report_url" in result and "report_URL" not in credentials:
                    credentials["report_URL"] = result["report_url"]

    # If export_URL not found in publishingResults, try to construct from main_report
    if "export_URL" not in credentials and "report_URL" in credentials:
        report_url = credentials["report_URL"]
        if "/3_Reports/analysis/" in report_url:
            credentials["export_URL"] = report_url.split("/3_Reports/analysis/")[0]

    return credentials
