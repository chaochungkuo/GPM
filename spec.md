# Add Export Engine API Integration to GPM

## Overview

Add support for sending export specifications to the export_engine API when the `--with-api` flag is used with the `gpm export` command. This converts GPM's export structure to the export_engine's `ExportJobSpec` format and submits it via HTTP POST.

## Implementation Steps

### 1. Add Configuration for Export Engine API

**File**: `config/gpm.ini`

- Add new `[EXPORT_ENGINE]` section with:
- `EXPORT_ENGINE_API_URL`: Base URL for the export engine API (e.g., `http://localhost:8000`)
- `EXPORT_ENGINE_BACKENDS`: Comma-separated list of backends to use (e.g., `apache,owncloud`)

### 2. Add CLI Flag to Export Command

**File**: `gpm/main.py`

- Add `--with-api` flag to the `export` command (around line 273)
- Pass this flag to the `pm.export()` method
- When `--with-api` is True, call a new method to send export spec to API

### 3. Create Export Specification Converter

**File**: `gpm/exports.py`

- Add function `convert_export_structure_to_job_spec()` that:
- Takes GPM's `export_structure` list, project profile, and export directory
- Converts each entry to `FileExport` format:
  - Maps source paths to absolute paths (using `self.prefix` for symlink handling)
  - Determines host the hostname command
  - Maps destination paths (target_dir + rename)
  - Sets mode to `symlink` (matching current GPM behavior)
  - Handles glob patterns by expanding them to individual FileExport entries
- Creates `ExportJobSpec` dict with:
  - `project_name`: from `self.profile["Project"]["project_name"]`
  - `export_list`: converted FileExport objects (as dicts)
  - `backend`: from config or default to `["apache"]`
  - `username`: from `self.profile["Export"]["export_user"]` if available
  - `password`: from `self.profile["Export"]["export_password"]` if available
  - `authors`: from `self.profile["Project"]["authors"]` if available (convert list to list of strings)

### 4. Add Host Detection Helper

**File**: `gpm/exports.py`

- Add function `determine_host_from_path()` that:
- Checks the hostname using the `hostname` command 
- Returns host name (e.g., `"nextgen"`, `"nextgen2"`, `"nextgen3"`)

### 5. Add API Submission Function

**File**: `gpm/exports.py`

- Add function `submit_export_to_api()` that:
- Takes the `ExportJobSpec` dict
- Reads API URL from config (`get_gpm_config("EXPORT_ENGINE", "EXPORT_ENGINE_API_URL")`)
- Makes POST request to `/export` endpoint using `requests.post()`
- Handles errors gracefully with user-friendly messages
- Returns job_id and status dict
- Uses `requests` library (already available via `helper.py`)

### 6. Integrate API Export into GPM Class

**File**: `gpm/gpm.py`

- Modify `export()` method to accept optional `use_api` parameter
- When `use_api=True`:
- Load export config (already done via `load_export_config()`)
- Convert export structure to job spec using `convert_export_structure_to_job_spec()`
- Submit to API using `submit_export_to_api()`
- Display job_id and status to user
- Optionally still perform local export (or make it mutually exclusive - need to clarify)
- Ensure source paths are resolved to absolute paths before conversion
- Handle symbolic link prefixes (already handled by `self.prefix`)

### 7. Error Handling and User Feedback

**Files**: `gpm/exports.py` and `gpm/gpm.py`

- Add try/except blocks around API calls
- Display clear error messages if API is unavailable or returns errors
- Show job_id and status when submission succeeds
- Provide helpful messages about checking job status

## Key Design Decisions

1. **Backward Compatibility**: The `--with-api` flag is optional, so existing workflows continue to work
2. **Path Handling**: Use existing `self.prefix` mechanism for symbolic link handling
3. **Host Detection**: First check for `/mnt/{host}` patterns, then fall back to hostname command, default to "nextgen"
4. **Mode**: Default to `symlink` mode to match current GPM behavior
5. **Backends**: Read from config, default to `["apache"]` if not specified
6. **API vs Local Export**: When `--with-api` is used, should we still do local export? (Need to clarify - for now, assume both can happen)

## Files to Modify

1. `config/gpm.ini` - Add export engine configuration section
2. `gpm/main.py` - Add `--with-api` CLI flag to export command
3. `gpm/exports.py` - Add API submission, conversion, and host detection functions
4. `gpm/gpm.py` - Integrate API export into export workflow

## Testing Considerations

- Test with valid export config
- Test with missing API URL (should fail gracefully)
- Test with invalid project paths
- Test host detection for different mount points (`/mnt/nextgen`, `/mnt/nextgen2`, `/data`)
- Test with and without `--with-api` flag
- Test with glob patterns in export config
- Test error handling for API failures