import glob
import json
import os
import subprocess
import sys
import threading
import time
from typing import Any, Dict, Optional, Tuple

import click
import requests
import websocket

from gpm.exports import generate_password
from gpm.helper import get_gpm_config


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


def determine_report_section(dest_path: str) -> str:
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
                            "src": matching_file,
                            "dest": dest_path,
                            "host": host,
                            "project": project_name,
                            "mode": "symlink",
                            "include_in_report": True,
                            "report_section": determine_report_section(dest_path),
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
                    "src": abs_source,
                    "dest": dest_path,
                    "host": host,
                    "project": project_name,
                    "include_in_report": True,
                    "mode": "symlink",
                    "report_section": determine_report_section(dest_path),
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
        if ws_url_config and isinstance(ws_url_config, str):
            ws_base_url = ws_url_config.rstrip("/")
        else:
            ws_base_url = _convert_api_url_to_ws_url(api_url)
    except Exception:
        ws_base_url = _convert_api_url_to_ws_url(api_url)

    ws_url = f"{ws_base_url}/export/ws/{job_id}"

    completion_data = None
    error_occurred = False
    premature_closure = False

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
            formatted_msg = notification.get("formatted_message", "")

            # Display messages to user
            if notif_type == "completion" and status == "done":
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
            else:
                # For other messages, display formatted_message if available, else message
                if formatted_msg:
                    click.echo(click.style(formatted_msg, fg="bright_blue"))
                elif msg:
                    click.echo(click.style(msg, fg="bright_blue"))

            # Handle completion message
            if notif_type == "completion" and status == "done":
                # Extract completion data from the notification
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
        nonlocal premature_closure
        if close_status_code and close_status_code != 1000:
            premature_closure = True
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
                premature_closure = True
        else:
            ws.run_forever()

        # If websocket closed prematurely and we don't have completion data, try polling
        if premature_closure and completion_data is None and not error_occurred:
            click.echo(
                click.style(
                    "\nWebSocket closed prematurely. Falling back to polling...",
                    fg="yellow",
                )
            )
            notification = poll_final_message(job_id, api_url)
            if notification:
                completion_data = {
                    "job_id": job_id,
                    "status": "done",
                    "notification": notification,
                }

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
        # If websocket connection fails, try polling as fallback
        click.echo(
            click.style(
                "WebSocket connection failed. Falling back to polling...",
                fg="yellow",
            )
        )
        notification = poll_final_message(job_id, api_url)
        if notification:
            return {
                "job_id": job_id,
                "status": "done",
                "notification": notification,
            }
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
