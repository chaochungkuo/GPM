"""
Unit tests for export_engine API integration in GPM.

Tests cover:
- Host determination
- Export structure to job spec conversion
- API submission
- Integration with GPM export workflow
"""

import os
import subprocess
import sys
from unittest.mock import Mock, patch

import pytest
import requests

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gpm.api_export import (
    convert_export_structure_to_job_spec,
    determine_host,
    submit_export_to_api,
)

# Static mock values used across multiple tests
DEFAULT_PASSWORD = "testpass123"
DEFAULT_HOST = "nextgen"
DEFAULT_BACKEND = "apache"
DEFAULT_API_URL = "http://localhost:8000"

PROFILE = {
    "Project": {
        "project_name": "240101_TestUser_Group_Institute_RNAseq",
        "project_path": "/mnt/nextgen/projects/test_project",
    }
}


@pytest.fixture
def mock_subprocess_run():
    """Fixture for mocking subprocess.run."""
    with patch("gpm.api_export.subprocess.run") as mock_run:
        yield mock_run


class TestDetermineHost:
    """Tests for determine_host() function."""

    def test_determine_host_success(self, mock_subprocess_run):
        """Test successful hostname determination."""
        mock_subprocess_run.return_value = Mock(stdout="nextgen\n", returncode=0)
        result = determine_host()
        assert result == "nextgen"

    def test_determine_host_with_domain(self, mock_subprocess_run):
        """Test hostname extraction when domain is present."""
        mock_subprocess_run.return_value = Mock(
            stdout="nextgen2.example.com\n", returncode=0
        )
        result = determine_host()
        assert result == "nextgen2"

    @patch("gpm.api_export.sys.exit")
    def test_determine_host_error(self, mock_exit, mock_subprocess_run):
        """Test error handling when hostname command fails."""
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(1, "hostname")
        determine_host()
        mock_exit.assert_called_once_with(1)


@pytest.fixture
def mock_export_structure(tmp_path):
    """Create a mock export structure with test files."""
    test_file = tmp_path / "test_file.txt"
    test_file.write_text("test content")
    return [
        ["RNAseq", str(test_file), "results", None],
    ]


@pytest.fixture
def mock_common_patches():
    """Fixture that patches common dependencies with default values."""
    with patch("gpm.exports.generate_password") as mock_password, patch(
        "gpm.api_export.get_gpm_config"
    ) as mock_config, patch("gpm.api_export.determine_host") as mock_host:
        mock_password.return_value = DEFAULT_PASSWORD
        mock_config.return_value = DEFAULT_BACKEND
        mock_host.return_value = DEFAULT_HOST
        yield {
            "password": mock_password,
            "config": mock_config,
            "host": mock_host,
        }


class TestConvertExportStructureToJobSpec:
    """Tests for convert_export_structure_to_job_spec() function."""

    def test_convert_basic_export_structure(
        self, mock_common_patches, mock_export_structure
    ):
        """Test basic conversion of export structure to job spec."""
        job_spec = convert_export_structure_to_job_spec(
            mock_export_structure, PROFILE, prefix=""
        )

        assert job_spec["project_name"] == PROFILE["Project"]["project_name"]
        assert (
            job_spec["username"]
            == PROFILE["Project"]["project_name"].split("_")[1].lower()
        )
        assert job_spec["password"] == DEFAULT_PASSWORD
        assert job_spec["backend"] == [DEFAULT_BACKEND]
        assert len(job_spec["export_list"]) == 1
        assert job_spec["export_list"][0]["mode"] == "symlink"
        assert job_spec["export_list"][0]["host"] == DEFAULT_HOST

    @patch("gpm.api_export.glob.glob")
    def test_convert_with_glob_pattern(self, mock_glob, mock_common_patches, tmp_path):
        """Test conversion with glob patterns in source paths."""
        # Create test files
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file1.write_text("content1")
        file2.write_text("content2")

        mock_glob.return_value = [str(file1), str(file2)]

        export_structure = [["RNAseq", "*.txt", "results", None]]
        job_spec = convert_export_structure_to_job_spec(
            export_structure, PROFILE, prefix=str(tmp_path) + "/"
        )

        assert len(job_spec["export_list"]) == 2
        assert all(entry["mode"] == "symlink" for entry in job_spec["export_list"])

    def test_convert_with_absolute_path(self, mock_common_patches, tmp_path):
        """Test conversion with absolute source paths."""
        test_file = tmp_path / "absolute_file.txt"
        test_file.write_text("content")

        export_structure = [["RNAseq", str(test_file), "results", "renamed.txt"]]
        job_spec = convert_export_structure_to_job_spec(
            export_structure, PROFILE, prefix=""
        )

        assert job_spec["export_list"][0]["source"] == str(test_file)
        assert job_spec["export_list"][0]["destination"] == "results/renamed.txt"

    def test_convert_with_authors(self, mock_common_patches, mock_export_structure):
        """Test conversion includes authors when available."""
        PROFILE["Project"]["authors"] = ["author1", "author2"]
        job_spec = convert_export_structure_to_job_spec(
            mock_export_structure, PROFILE, prefix=""
        )

        assert "authors" in job_spec
        assert job_spec["authors"] == ["author1", "author2"]

    def test_convert_backend_from_config(
        self, mock_common_patches, mock_export_structure
    ):
        """Test backend configuration from config file."""
        mock_common_patches["config"].return_value = "apache,owncloud"

        job_spec = convert_export_structure_to_job_spec(
            mock_export_structure, PROFILE, prefix=""
        )

        assert job_spec["backend"] == ["apache", "owncloud"]

    def test_convert_username_extraction(
        self, mock_common_patches, mock_export_structure
    ):
        """Test username extraction from project name."""
        # Test with standard format
        job_spec = convert_export_structure_to_job_spec(
            mock_export_structure, PROFILE, prefix=""
        )
        assert job_spec["username"] == "testuser"

        # Test with non-standard format (fallback)
        PROFILE["Project"]["project_name"] = "InvalidName"
        job_spec = convert_export_structure_to_job_spec(
            mock_export_structure, PROFILE, prefix=""
        )
        assert job_spec["username"] == "invalidname"


@pytest.fixture
def mock_job_spec():
    """Create a mock job specification."""
    return {
        "project_name": "test_project",
        "export_list": [{"source": "/path/to/file", "destination": "results/file"}],
        "backend": ["apache"],
        "username": "testuser",
        "password": "testpass",
    }


@pytest.fixture
def mock_api_patches():
    """Fixture that patches API-related dependencies."""
    with patch("gpm.api_export.get_gpm_config") as mock_config, patch(
        "gpm.api_export.requests.post"
    ) as mock_post:
        mock_config.return_value = DEFAULT_API_URL
        yield {
            "config": mock_config,
            "post": mock_post,
        }


class TestSubmitExportToAPI:
    """Tests for submit_export_to_api() function."""

    def test_submit_success(self, mock_api_patches, mock_job_spec):
        """Test successful API submission."""
        mock_response = Mock()
        mock_response.status_code = 201
        mock_response.json.return_value = {"job_id": "job123", "status": "submitted"}
        mock_api_patches["post"].return_value = mock_response

        job_id, result = submit_export_to_api(mock_job_spec)

        assert job_id == "job123"
        assert isinstance(result, dict)
        assert result.get("status") == "submitted"
        mock_api_patches["post"].assert_called_once()
        call_args = mock_api_patches["post"].call_args
        assert call_args[1]["json"] == mock_job_spec

    def test_submit_api_error(self, mock_api_patches, mock_job_spec):
        """Test API error response handling."""
        mock_response = Mock()
        mock_response.status_code = 400
        mock_response.json.return_value = {"error": "Invalid request"}
        mock_response.text = '{"error": "Invalid request"}'
        mock_api_patches["post"].return_value = mock_response

        job_id, error = submit_export_to_api(mock_job_spec)

        assert job_id is None
        assert "400" in error
        assert "Invalid request" in error

    def test_submit_connection_error(self, mock_api_patches, mock_job_spec):
        """Test connection error handling."""
        mock_api_patches["post"].side_effect = requests.exceptions.ConnectionError(
            "Connection refused"
        )

        job_id, error = submit_export_to_api(mock_job_spec)

        assert job_id is None
        assert "Failed to connect" in error

    @patch("gpm.api_export.get_gpm_config")
    def test_submit_missing_config(self, mock_config, mock_job_spec):
        """Test handling when API URL is not configured."""
        mock_config.return_value = None

        job_id, error = submit_export_to_api(mock_job_spec)

        assert job_id is None
        assert "not configured" in error


@pytest.fixture
def mock_gpm():
    """Create a mock GPM instance."""
    from gpm.gpm import GPM

    gpm = GPM()
    # Initialize profile as OrderedDict to match GPM structure
    gpm.profile["Project"]["project_name"] = "240101_TestUser_Group_Institute_RNAseq"
    gpm.profile["Project"]["project_path"] = "/mnt/nextgen/projects/test"
    gpm.export_structure = [
        ["RNAseq", "/path/to/file.txt", "results", None],
    ]
    gpm.prefix = ""
    gpm.load_export_config = Mock()
    gpm.update_project_name = Mock()
    return gpm


@pytest.fixture
def mock_export_patches():
    """Fixture that patches export-related dependencies."""
    with patch("gpm.gpm.check_export_directory") as mock_check_dir, patch(
        "gpm.gpm.convert_export_structure_to_job_spec"
    ) as mock_convert, patch("gpm.gpm.submit_export_to_api") as mock_submit, patch(
        "gpm.gpm.click.echo"
    ) as mock_echo:
        mock_check_dir.return_value = None
        yield {
            "check_dir": mock_check_dir,
            "convert": mock_convert,
            "submit": mock_submit,
            "echo": mock_echo,
        }


class TestExportIntegration:
    """Integration tests for export() method with API support."""

    def test_export_with_api_success(self, mock_export_patches, mock_gpm, tmp_path):
        """Test export with API flag when submission succeeds."""
        export_dir = str(tmp_path / "export")
        mock_export_patches["convert"].return_value = {"project_name": "test"}
        mock_export_patches["submit"].return_value = ("job123", {"status": "submitted"})

        mock_gpm.export(export_dir, use_api=True, symlink=False)

        mock_export_patches["convert"].assert_called_once()
        mock_export_patches["submit"].assert_called_once()
        # Verify success message was printed
        assert any(
            "successfully" in str(call).lower()
            for call in mock_export_patches["echo"].call_args_list
        )

    def test_export_with_api_error(self, mock_export_patches, mock_gpm, tmp_path):
        """Test export with API flag when submission fails."""
        export_dir = str(tmp_path / "export")
        mock_export_patches["convert"].return_value = {"project_name": "test"}
        mock_export_patches["submit"].return_value = (None, "API error message")

        mock_gpm.export(export_dir, use_api=True, symlink=False)

        mock_export_patches["submit"].assert_called_once()
        # Verify error message was printed
        assert any(
            "failed" in str(call).lower()
            for call in mock_export_patches["echo"].call_args_list
        )
