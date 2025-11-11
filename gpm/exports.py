import glob
import hashlib
import json
import os
import random
import re
import string
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import click
import requests
import websocket
import xtarfile as tarfile
from owncloud import Client  # type: ignore
from tqdm import tqdm

from gpm.helper import get_gpm_config, get_gpmdata_path


def check_export_directory(export_folder):
    # Check if the folder exists
    if not os.path.exists(export_folder):
        # If not, create the folder
        os.makedirs(export_folder)
        print(f"Folder '{export_folder}' created.")
    else:
        print(f"Folder '{export_folder}' already exists.")
        sys.exit()


def get_htaccess_path():
    """
    Return the path of htaccess
    """
    config_path = os.path.join(get_gpmdata_path(), "config/htaccess")
    return config_path


def htpasswd_create_user(export_folder, url, username, app=None):
    """Create the new user in the target directory with password"""
    export_base_path = Path(export_folder).parent.absolute()
    if os.path.exists(os.path.join(export_base_path, ".htpasswd")):
        subprocess.run(
            [
                "cp",
                os.path.join(export_base_path, ".htpasswd"),
                os.path.join(export_folder, ".htpasswd"),
            ]
        )
        # shutil.copy(os.path.join(export_base_path, ".htpasswd"),
        #             os.path.join(export_folder, ".htpasswd"))
        password = generate_password()
        cmd = " ".join(
            [
                "htpasswd",
                "-b",
                os.path.join(export_folder, ".htpasswd"),
                username,
                password,
            ]
        )
        subprocess.run(cmd, shell=True)
        click.echo()
        click.echo(
            click.style("Create new user for export directory:", fg="bright_green")
        )
        click.echo("Directory:\t" + export_folder)
        export_URL = url
        if app:
            if app in ["RNAseq", "tRNAseq", "mRNAseq", "3mRNAseq"]:
                app = "RNAseq"
            for repo_app in get_gpm_config("GPM", "GPM_REPORTS"):
                if repo_app.lower() == app.lower():
                    export_URL = "".join(
                        [url, "/3_Reports/analysis/Analysis_Report_", repo_app, ".html"]
                    )
        # click.echo("user:\t" + username)
        # click.echo("password:\t" + password)
        # click.echo("URL:\t" + export_URL)
        return username, password
    else:
        click.echo("Skip setting htpasswd")
        return None, None


def create_user(export_folder, export_URL, username):
    htpasswd_create_user(export_folder, export_URL, username, None)


def generate_password():
    source = string.ascii_letters + string.digits
    result_str = "".join((random.choice(source) for i in range(12)))
    return result_str


def export_empty_folder(export_URL, export_dir, username):
    # Add htaccess
    # data_dir = os.path.join(os.path.dirname(__file__), "data")
    htaccess_path = get_htaccess_path()
    with open(htaccess_path) as f1:
        contents = [le.rstrip() for le in f1.readlines()]
    for i, line in enumerate(contents):
        if "GPM_TITLE_NAME" in line:
            contents[i] = line.replace("GPM_TITLE_NAME", os.path.basename(export_dir))

    with open(os.path.join(export_dir, ".htaccess"), "w") as f2:
        for line in contents:
            if "PROJECT_PROJECT_NAME" in line:
                line = line.replace(
                    "PROJECT_PROJECT_NAME", os.path.basename(export_dir)
                )
            print(line, file=f2)
    # Create user
    _, password = htpasswd_create_user(
        export_dir,
        os.path.join(export_URL, os.path.basename(export_dir)),
        username.lower(),
        None,
    )

    oc = owncloud_login()
    url = owncloud_export(oc, os.path.basename(export_dir), password)
    # click.echo("Download URL: " + url)


def owncloud_login() -> Client:
    owncloud_pass = os.getenv("OWNCLOUD_SHARE_PASS")
    if owncloud_pass is None:
        click.echo(
            "Could not find the password for the owncloud account, please set the environment variable `OWNCLOUD_SHARE_PASS`"
        )
        sys.exit()
    oc = Client("https://genomics.rwth-aachen.de/cloud")
    oc.login("GPM", owncloud_pass)
    return oc


def owncloud_export(oc, path, password=None) -> str:
    PATH_PREFIX = get_gpm_config("EXPORT", "EXPORT_CLOUD_PREFIX")
    out = oc.share_file_with_link(
        str(Path(PATH_PREFIX, path)), password=password, public_upload=False
    )
    return out.get_link()


def relpath(path_file):
    base_dirs = ["/mnt/nextgen", "/mnt/nextgen2", "/mnt/nextgen3"]
    for base_dir in base_dirs:
        if path_file.startswith(base_dir):
            # Getting the relative path of the directory
            rel_path = os.path.relpath(path_file, base_dir)
            print("rel_path: " + rel_path)
            path_file = os.path.join("/", rel_path)
            print("path_file: " + path_file)
    return path_file


def tar_exports(export_folder, dry_run, gzip, same_server=False):
    def fits_pattern(filename, pattern):
        return bool(re.match(pattern.replace("*", ".*"), filename))

    regex_patterns = get_gpm_config("EXPORT", "TAR_EXPORT_IGNORE")
    if export_folder == ".":
        export_folder = os.getcwd()
    export_folder = export_folder.rstrip("/")
    name = os.path.basename(export_folder)
    compressed_folder = os.path.join(export_folder, "compressed_tar")
    # Create compressed_tar folder
    if not os.path.exists(compressed_folder):
        click.echo(click.style("Create the folder:", fg="bright_green"))
        click.echo(compressed_folder)
        if not dry_run:
            os.makedirs(compressed_folder)
    # Tar each folder
    for filename in os.listdir(export_folder):
        if filename.startswith("."):
            continue
        fits_any_pattern = any(
            fits_pattern(filename, pattern) for pattern in regex_patterns
        )
        if fits_any_pattern:
            continue
        else:
            path_file = os.path.join(export_folder, filename)
            if gzip:
                tarfile = os.path.join(
                    compressed_folder, name + "_" + filename + ".tar.gz"
                )
            else:
                tarfile = os.path.join(
                    compressed_folder, name + "_" + filename + ".tar"
                )
            # print("path_file: " + path_file)
            if os.path.islink(path_file):
                path_file = os.readlink(path_file)
                # print("path_file link: " + path_file)
                # path_file = relpath(path_file)
                print("path_file link: " + path_file)
            if os.path.isdir(path_file) and "compressed_tar" not in filename:
                click.echo(click.style("Tar the folder:", fg="bright_green"))
                click.echo(path_file + click.style(" => ", fg="bright_green") + tarfile)
                if not dry_run:
                    tar_folder(path_file, tarfile, gzip)


def tar_folder(input_folder, output_tar, gzip):
    regex_patterns = get_gpm_config("EXPORT", "TAR_EXPORT_IGNORE")
    if gzip:
        tar_mode = "w:gz"
    else:
        tar_mode = "w"
    if os.path.exists(output_tar):
        click.echo(output_tar + " exists.")
    else:
        with tarfile.open(output_tar, tar_mode) as tar:
            # Set up the tqdm progress bar
            total_size = 0
            for root, dirs, files in os.walk(input_folder, followlinks=True):
                # Filter root
                if os.path.basename(root) in regex_patterns:
                    continue
                # Filter directories
                for pattern in regex_patterns:
                    if pattern in dirs:
                        dirs.remove(pattern)
                        continue
                # Iterate files
                for name in files:
                    full_path = os.path.join(root, name)
                    total_size += get_size(full_path)
            progress_bar = tqdm(
                total=total_size, desc="Creating tar archive", unit="B", unit_scale=True
            )
            for root, dirs, files in os.walk(input_folder, followlinks=True):
                # Filter root
                if os.path.basename(root) in regex_patterns:
                    continue
                for pattern in regex_patterns:
                    if pattern in dirs:
                        dirs.remove(pattern)  # Ignore 'renv' folder
                        continue
                for name in files:
                    full_path = os.path.join(root, name)
                    arcname = os.path.relpath(full_path, input_folder)
                    # Add the file or directory to the tar archive
                    tar.add(full_path, arcname=arcname)
                    # Update the progress bar by the size
                    progress_bar.update(get_size(full_path))
            # Close the progress bar
            progress_bar.close()
        save_md5_to_file(output_tar)


def get_size(path):
    # If path is a symbolic link, get the size of the target file
    if os.path.islink(path):
        target_path = os.path.realpath(path)
        target_stat = os.stat(target_path)
        return target_stat.st_size
    elif os.path.isdir(path):
        return get_folder_size(path)
    else:
        return os.path.getsize(path)


def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            md5.update(chunk)
    return md5.hexdigest()


def save_md5_to_file(file_path):
    md5_checksum = calculate_md5(file_path)
    with open(file_path + ".md5", "w") as md5_file:
        md5_file.write(md5_checksum)


def get_folder_size(folder_path):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(folder_path):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            total_size += os.path.getsize(file_path)
    return total_size


def determine_host():
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


def convert_export_structure_to_job_spec(export_structure, profile, prefix=""):
    """
    Convert GPM's export structure to export_engine's ExportJobSpec format.

    Args:
        export_structure: List of export entries, each with [app, source, target_dir, rename]
        profile: Project profile dictionary
        export_dir: Target export directory
        prefix: Symbolic link prefix (for handling symlinks)

    Returns:
        dict: ExportJobSpec in the format expected by the export engine API
    """
    project_name = profile["Project"]["project_name"]
    project_path = profile["Project"]["project_path"]

    # Generate username from project_name (format: date_username_group_institute_application)
    # Extract the second component (index 1) as username
    project_parts = project_name.split("_")
    if len(project_parts) >= 2:
        username = project_parts[1].lower()
    else:
        # Fallback if project_name doesn't match expected format
        username = project_name.lower()

    # Generate password using existing function
    password = generate_password()

    # Get backend from config or default to ["apache"]
    try:
        backends_str = get_gpm_config("EXPORT_ENGINE", "EXPORT_ENGINE_BACKENDS")
        if isinstance(backends_str, str):
            backend = [b.strip() for b in backends_str.split(",")]
        else:
            backend = backends_str if isinstance(backends_str, list) else ["apache"]
    except Exception:
        backend = ["apache"]

    # Get authors if available
    authors = []
    if "authors" in profile["Project"] and profile["Project"]["authors"]:
        authors_list = profile["Project"]["authors"]
        if isinstance(authors_list, list):
            authors = [str(a) for a in authors_list]
        elif isinstance(authors_list, str):
            authors = [authors_list]

    # Determine host
    host = determine_host()

    # Convert export structure to FileExport format
    export_list = []
    for entry in export_structure:
        # entry format: [app, source_path, target_dir, rename]
        if not entry[1]:  # Skip entries without source (folder creation only)
            continue

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
            # Expand glob pattern
            if os.path.isabs(source_path):
                glob_pattern = prefix + source_path
            else:
                glob_pattern = prefix + os.path.join(project_path, source_path)

            matching_files = glob.glob(glob_pattern)
            for matching_file in matching_files:
                # Only include files (not directories) from glob patterns
                # This matches the original export behavior
                if os.path.isfile(matching_file):
                    # Determine destination filename
                    if rename:
                        dest_filename = rename
                    else:
                        dest_filename = os.path.basename(matching_file)

                    dest_path = os.path.join(target_dir, dest_filename)

                    export_list.append(
                        {
                            "source": matching_file,
                            "destination": dest_path,
                            "host": host,
                            "mode": "symlink",
                        }
                    )
        else:
            # Single file or directory
            # Include in spec even if it doesn't exist yet - export engine will handle
            # Determine destination path
            if rename:
                dest_path = os.path.join(target_dir, rename)
            else:
                dest_path = os.path.join(target_dir, os.path.basename(abs_source))

            export_list.append(
                {
                    "source": abs_source,
                    "destination": dest_path,
                    "host": host,
                    "mode": "symlink",
                }
            )

    # Create ExportJobSpec dict
    job_spec = {
        "project_name": project_name,
        "export_list": export_list,
        "backend": backend,
        "username": username,
        "password": password,
    }

    # Add authors if available
    if authors:
        job_spec["authors"] = authors

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
        api_url = api_url.rstrip("/")
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


def _convert_api_url_to_ws_url(api_url: str) -> str:
    """
    Convert HTTP/HTTPS API URL to WebSocket URL (WS/WSS).

    Args:
        api_url: Base API URL (e.g., "http://localhost:8000" or "https://api.example.com")

    Returns:
        WebSocket URL (e.g., "ws://localhost:8000" or "wss://api.example.com")
    """
    api_url = api_url.rstrip("/")
    if api_url.startswith("https://"):
        return api_url.replace("https://", "wss://", 1)
    elif api_url.startswith("http://"):
        return api_url.replace("http://", "ws://", 1)
    else:
        # Assume HTTP if no scheme
        return f"ws://{api_url}"


def monitor_job_via_websocket(
    job_id: str, api_url: str, timeout: Optional[float] = None
) -> Optional[Dict[str, Any]]:
    """
    Monitor export job via WebSocket for real-time status updates.

    Args:
        job_id: The job ID to monitor
        api_url: Base API URL (will be converted to WebSocket URL)
        timeout: Optional timeout in seconds (None = no timeout)

    Returns:
        Dictionary with completion data (credentials, URLs) if job completes successfully,
        None if job fails or connection error occurs
    """
    # Get WebSocket URL from config or convert from API URL
    try:
        ws_url_config = get_gpm_config("EXPORT_ENGINE", "EXPORT_ENGINE_WS_URL")
        if ws_url_config:
            ws_base_url = ws_url_config.rstrip("/")
        else:
            ws_base_url = _convert_api_url_to_ws_url(api_url)
    except Exception:
        ws_base_url = _convert_api_url_to_ws_url(api_url)

    ws_url = f"{ws_base_url}/export/ws/{job_id}"

    completion_data = None
    error_occurred = False

    def on_message(ws, message):
        """Handle incoming WebSocket messages."""
        nonlocal completion_data, error_occurred

        try:
            # Handle ping/pong keepalive - server sends "ping" as text, we respond with "pong"
            if message == "ping":
                ws.send("pong")
                return
            elif message == "pong":
                return

            # Parse JSON notification
            try:
                notification = json.loads(message)
            except json.JSONDecodeError:
                click.echo(
                    click.style(
                        f"Warning: Received invalid JSON message: {message[:100]}",
                        fg="yellow",
                    )
                )
                return

            # Extract notification fields
            notif_type = notification.get("type", "normal")
            status = notification.get("status", "")
            msg = notification.get("message", "")
            formatted_msg = notification.get("formatted_message", msg)

            # Display formatted message to user
            if formatted_msg:
                click.echo(click.style(formatted_msg, fg="bright_blue"))
            elif msg:
                click.echo(click.style(msg, fg="bright_blue"))

            # Handle completion message
            if notif_type == "completion" and status == "done":
                # Extract completion data from the formatted message or notification
                # The server sends publisher results in the completion message
                completion_data = {
                    "job_id": job_id,
                    "status": status,
                    "notification": notification,
                }
                # Close WebSocket after receiving completion
                ws.close(status=1000, reason="Job completed successfully")

            # Handle error message
            elif notif_type == "error" or status == "failed":
                error_occurred = True
                click.echo(click.style(f"Job failed: {msg}", fg="red"))
                ws.close(status=1000, reason="Job failed")

        except Exception as e:
            click.echo(
                click.style(f"Error processing WebSocket message: {str(e)}", fg="red")
            )
            error_occurred = True
            ws.close(status=1000, reason="Error processing message")

    def on_error(ws, error):
        """Handle WebSocket errors."""
        nonlocal error_occurred
        error_occurred = True
        click.echo(click.style(f"WebSocket error: {str(error)}", fg="red"))

    def on_close(ws, close_status_code, close_msg):
        """Handle WebSocket close."""
        if close_status_code and close_status_code != 1000:
            click.echo(
                click.style(
                    f"WebSocket closed unexpectedly: {close_status_code}", fg="yellow"
                )
            )

    def on_open(ws):
        """Handle WebSocket open."""
        click.echo(
            click.style(
                "Connected to export engine. Monitoring job progress...",
                fg="bright_green",
            )
        )

    # Create WebSocket connection
    ws = None
    try:
        ws = websocket.WebSocketApp(
            ws_url,
            on_message=on_message,
            on_error=on_error,
            on_close=on_close,
            on_open=on_open,
        )

        # Run WebSocket with optional timeout
        # Note: websocket-client's run_forever() doesn't support timeout directly,
        # so we use a thread-based approach or run in a separate thread with timeout
        import threading

        if timeout:

            def run_with_timeout():
                ws.run_forever()

            thread = threading.Thread(target=run_with_timeout, daemon=True)
            thread.start()
            thread.join(timeout=timeout)

            if thread.is_alive():
                click.echo(
                    click.style(
                        f"Monitoring timeout after {timeout} seconds", fg="yellow"
                    )
                )
                ws.close(status=1000, reason="Monitoring timeout")
                return None
        else:
            ws.run_forever()

        return completion_data if not error_occurred else None

    except KeyboardInterrupt:
        click.echo(click.style("\nMonitoring interrupted by user", fg="yellow"))
        if ws is not None:
            try:
                ws.close(status=1000, reason="Monitoring interrupted by user")
            except Exception:
                pass  # Ignore errors when closing
        return None
    except Exception as e:
        click.echo(click.style(f"Failed to connect to WebSocket: {str(e)}", fg="red"))
        if ws is not None:
            try:
                ws.close(status=1000, reason="Connection error")
            except Exception:
                pass  # Ignore errors when closing
        return None


def extract_credentials_from_completion(
    completion_data: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Extract credentials and URLs from completion notification message.

    The completion message contains a formatted_message that includes credentials.
    We parse the plain message which has a structured format.

    Args:
        completion_data: Dictionary containing completion notification

    Returns:
        Dictionary with extracted credentials: username, password, export_URL, report_URL, download_url
    """
    if not completion_data or "notification" not in completion_data:
        return {}

    notification = completion_data["notification"]
    message = notification.get("message", "")

    # The plain message has a structured format like:
    # 'Project ID': '...', 'Report URL': '...', 'Username': '...', 'Password': '...', etc.
    credentials = {}

    # Extract from plain message using regex
    username_match = re.search(r"'Username':\s*'([^']+)'", message)
    password_match = re.search(r"'Password':\s*'([^']+)'", message)
    report_url_match = re.search(r"'Report URL':\s*'([^']+)'", message)
    download_url_match = re.search(r"'Download URL':\s*'([^']+)'", message)

    if username_match:
        credentials["username"] = username_match.group(1)
    if password_match:
        credentials["password"] = password_match.group(1)
    if report_url_match:
        credentials["report_URL"] = report_url_match.group(1)
    if download_url_match:
        credentials["download_url"] = download_url_match.group(1)

    # Also try to extract export URL from formatted message or publisher results
    # The formatted message may contain publisher URLs
    formatted_msg = notification.get("formatted_message", "")

    # Look for Apache URL in formatted message
    apache_url_match = re.search(r"\*\*URL:\*\*\s*([^\s\n]+)", formatted_msg)
    if apache_url_match:
        credentials["export_URL"] = apache_url_match.group(1)
    elif "export_URL" not in credentials:
        # Fallback: construct from report URL if available
        if "report_URL" in credentials:
            # Extract base URL from report URL
            report_url = credentials["report_URL"]
            if "/3_Reports/analysis/" in report_url:
                credentials["export_URL"] = report_url.split("/3_Reports/analysis/")[0]

    return credentials
