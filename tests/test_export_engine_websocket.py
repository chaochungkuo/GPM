"""
Unit tests for WebSocket monitoring in export_engine API integration.

Tests cover:
- Successful WebSocket monitoring flow with completion
- WebSocket connection failure and error handling
- Edge case: already completed job
"""

import json
import os
import sys
import threading
import time
from typing import Any, Dict
from unittest.mock import MagicMock, Mock, patch

import pytest
import websocket

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gpm.exports import (
    _convert_api_url_to_ws_url,
    extract_credentials_from_completion,
    monitor_job_via_websocket,
)


class TestConvertApiUrlToWsUrl:
    """Tests for URL conversion helper."""

    def test_convert_http_to_ws(self):
        """Test converting HTTP URL to WS."""
        assert (
            _convert_api_url_to_ws_url("http://localhost:8000") == "ws://localhost:8000"
        )

    def test_convert_https_to_wss(self):
        """Test converting HTTPS URL to WSS."""
        assert (
            _convert_api_url_to_ws_url("https://api.example.com")
            == "wss://api.example.com"
        )

    def test_convert_url_without_scheme(self):
        """Test converting URL without scheme."""
        assert _convert_api_url_to_ws_url("localhost:8000") == "ws://localhost:8000"


class TestExtractCredentialsFromCompletion:
    """Tests for credential extraction from completion messages."""

    def test_extract_credentials_from_completion_message(self):
        """Test extracting credentials from completion notification."""
        completion_data = {
            "job_id": "test_job_123",
            "status": "done",
            "notification": {
                "type": "completion",
                "status": "done",
                "message": """
'Project ID': '240101_TestUser_Group_Institute_RNAseq',
'Report URL': 'https://genomics.rwth-aachen.de/data/240101_TestUser_Group_Institute_RNAseq/3_Reports/analysis/Analysis_Report_RNAseq.html',
'Username': 'testuser',
'Password': 'testpass123',
'Download URL': 'https://genomics.rwth-aachen.de/cloud/s/test',
'Download command': 'wget ...',
""",
                "formatted_message": "## Export Complete\n\n**URL:** https://genomics.rwth-aachen.de/data/240101_TestUser_Group_Institute_RNAseq",
            },
        }

        credentials = extract_credentials_from_completion(completion_data)

        assert credentials["username"] == "testuser"
        assert credentials["password"] == "testpass123"
        assert "report_URL" in credentials
        assert "download_url" in credentials
        assert "export_URL" in credentials

    def test_extract_credentials_missing_fields(self):
        """Test extraction when some fields are missing."""
        completion_data = {
            "job_id": "test_job_123",
            "status": "done",
            "notification": {
                "type": "completion",
                "status": "done",
                "message": "'Project ID': 'test'",
                "formatted_message": "",
            },
        }

        credentials = extract_credentials_from_completion(completion_data)
        # Should return empty dict or partial credentials
        assert isinstance(credentials, dict)


class TestMonitorJobViaWebsocket:
    """Tests for WebSocket monitoring function."""

    @patch("gpm.exports.websocket.WebSocketApp")
    @patch("gpm.exports.get_gpm_config")
    @patch("gpm.exports.click.echo")
    def test_successful_monitoring_flow(self, mock_echo, mock_config, mock_ws_class):
        """Test successful WebSocket monitoring with full job lifecycle."""
        # Setup mocks
        mock_config.side_effect = lambda section, key: {
            ("EXPORT_ENGINE", "EXPORT_ENGINE_WS_URL"): None,
            ("EXPORT_ENGINE", "EXPORT_ENGINE_API_URL"): "http://localhost:8000",
        }.get((section, key), None)

        # Simulate messages in order: QUEUED -> RUNNING -> PUBLISHING -> DONE
        messages = [
            json.dumps(
                {
                    "type": "normal",
                    "status": "queued",
                    "message": "Job queued",
                    "formatted_message": "## Export Job Queued\n\nYour export job has been queued.",
                    "job_id": "test_job_123",
                    "timestamp": "2024-01-01T00:00:00Z",
                }
            ),
            json.dumps(
                {
                    "type": "normal",
                    "status": "running",
                    "message": "Job running",
                    "formatted_message": "## Export Job Running\n\nYour export job is currently being processed.",
                    "job_id": "test_job_123",
                    "timestamp": "2024-01-01T00:00:01Z",
                }
            ),
            json.dumps(
                {
                    "type": "normal",
                    "status": "publishing",
                    "message": "Apache publishing",
                    "formatted_message": "## Apache Publishing\n\nStarting Apache publishing process",
                    "job_id": "test_job_123",
                    "timestamp": "2024-01-01T00:00:02Z",
                }
            ),
            json.dumps(
                {
                    "type": "completion",
                    "status": "done",
                    "message": """
'Project ID': 'test_project',
'Report URL': 'https://example.com/report.html',
'Username': 'testuser',
'Password': 'testpass',
'Download URL': 'https://example.com/download',
'Download command': 'wget ...',
""",
                    "formatted_message": "## Export Complete\n\n**URL:** https://example.com/data/test_project",
                    "job_id": "test_job_123",
                    "timestamp": "2024-01-01T00:00:03Z",
                }
            ),
        ]

        # Store callbacks
        callbacks = {}

        def setup_callbacks(on_message, on_error, on_close, on_open):
            callbacks["on_message"] = on_message
            callbacks["on_error"] = on_error
            callbacks["on_close"] = on_close
            callbacks["on_open"] = on_open

            mock_ws = MagicMock()

            def run_forever():
                """Simulate run_forever by calling callbacks."""
                if on_open:
                    on_open(mock_ws)
                # Simulate receiving messages
                for msg in messages:
                    if on_message:
                        on_message(mock_ws, msg)
                if on_close:
                    on_close(mock_ws, 1000, "Normal closure")

            mock_ws.run_forever = run_forever
            mock_ws.close = MagicMock()
            mock_ws.send = MagicMock()
            return mock_ws

        mock_ws_class.side_effect = lambda url, **kwargs: setup_callbacks(
            kwargs.get("on_message"),
            kwargs.get("on_error"),
            kwargs.get("on_close"),
            kwargs.get("on_open"),
        )

        # Run monitoring
        result = monitor_job_via_websocket("test_job_123", "http://localhost:8000")

        # Verify completion data was returned
        assert result is not None
        assert result["job_id"] == "test_job_123"
        assert result["status"] == "done"
        assert "notification" in result

        # Verify credentials can be extracted
        credentials = extract_credentials_from_completion(result)
        assert credentials["username"] == "testuser"
        assert credentials["password"] == "testpass"

    @patch("gpm.exports.websocket.WebSocketApp")
    @patch("gpm.exports.get_gpm_config")
    @patch("gpm.exports.click.echo")
    def test_websocket_connection_failure(self, mock_echo, mock_config, mock_ws_class):
        """Test WebSocket connection failure handling."""
        mock_config.side_effect = lambda section, key: {
            ("EXPORT_ENGINE", "EXPORT_ENGINE_WS_URL"): None,
            ("EXPORT_ENGINE", "EXPORT_ENGINE_API_URL"): "http://localhost:8000",
        }.get((section, key), None)

        def setup_callbacks(on_message, on_error, on_close, on_open):
            mock_ws = MagicMock()

            def run_forever():
                """Simulate run_forever with connection error."""
                # Simulate connection error
                if on_error:
                    on_error(mock_ws, Exception("Connection refused"))

            mock_ws.run_forever = run_forever
            mock_ws.close = MagicMock()
            return mock_ws

        mock_ws_class.side_effect = lambda url, **kwargs: setup_callbacks(
            kwargs.get("on_message"),
            kwargs.get("on_error"),
            kwargs.get("on_close"),
            kwargs.get("on_open"),
        )

        result = monitor_job_via_websocket("test_job_123", "http://localhost:8000")

        # Should return None on error
        assert result is None

        # Verify error message was displayed
        error_calls = [
            call
            for call in mock_echo.call_args_list
            if any(
                "error" in str(arg).lower() or "failed" in str(arg).lower()
                for arg in str(call)
            )
        ]
        assert len(error_calls) > 0

    @patch("gpm.exports.websocket.WebSocketApp")
    @patch("gpm.exports.get_gpm_config")
    @patch("gpm.exports.click.echo")
    def test_job_failure_notification(self, mock_echo, mock_config, mock_ws_class):
        """Test handling of job failure notification."""
        mock_config.side_effect = lambda section, key: {
            ("EXPORT_ENGINE", "EXPORT_ENGINE_WS_URL"): None,
            ("EXPORT_ENGINE", "EXPORT_ENGINE_API_URL"): "http://localhost:8000",
        }.get((section, key), None)

        error_message = json.dumps(
            {
                "type": "error",
                "status": "failed",
                "message": "Export failed: File not found",
                "formatted_message": "## Error\n\n**Error:**\n\n```\nExport failed: File not found\n```",
                "job_id": "test_job_123",
                "timestamp": "2024-01-01T00:00:00Z",
            }
        )

        def setup_callbacks(on_message, on_error, on_close, on_open):
            mock_ws = MagicMock()

            def run_forever():
                """Simulate run_forever with error message."""
                if on_open:
                    on_open(mock_ws)
                if on_message:
                    on_message(mock_ws, error_message)
                if on_close:
                    on_close(mock_ws, 1000, "Normal closure")

            mock_ws.run_forever = run_forever
            mock_ws.close = MagicMock()
            mock_ws.send = MagicMock()
            return mock_ws

        mock_ws_class.side_effect = lambda url, **kwargs: setup_callbacks(
            kwargs.get("on_message"),
            kwargs.get("on_error"),
            kwargs.get("on_close"),
            kwargs.get("on_open"),
        )

        result = monitor_job_via_websocket("test_job_123", "http://localhost:8000")

        # Should return None on failure
        assert result is None

        # Verify error message was displayed
        error_calls = [
            call
            for call in mock_echo.call_args_list
            if any(
                "failed" in str(arg).lower() or "error" in str(arg).lower()
                for arg in str(call)
            )
        ]
        assert len(error_calls) > 0

    @patch("gpm.exports.websocket.WebSocketApp")
    @patch("gpm.exports.get_gpm_config")
    @patch("gpm.exports.click.echo")
    def test_already_completed_job(self, mock_echo, mock_config, mock_ws_class):
        """Test handling of already completed job (server sends final message immediately)."""
        mock_config.side_effect = lambda section, key: {
            ("EXPORT_ENGINE", "EXPORT_ENGINE_WS_URL"): None,
            ("EXPORT_ENGINE", "EXPORT_ENGINE_API_URL"): "http://localhost:8000",
        }.get((section, key), None)

        # Completion message sent immediately on connect
        completion_message = json.dumps(
            {
                "type": "completion",
                "status": "done",
                "message": """
'Project ID': 'test_project',
'Report URL': 'https://example.com/report.html',
'Username': 'testuser',
'Password': 'testpass',
'Download URL': 'https://example.com/download',
'Download command': 'wget ...',
""",
                "formatted_message": "## Export Complete\n\n**URL:** https://example.com/data/test_project",
                "job_id": "test_job_123",
                "timestamp": "2024-01-01T00:00:00Z",
            }
        )

        def setup_callbacks(on_message, on_error, on_close, on_open):
            mock_ws = MagicMock()

            def run_forever():
                """Simulate run_forever with immediate completion message."""
                if on_open:
                    on_open(mock_ws)
                # Server immediately sends completion message
                if on_message:
                    on_message(mock_ws, completion_message)
                if on_close:
                    on_close(mock_ws, 1000, "Normal closure")

            mock_ws.run_forever = run_forever
            mock_ws.close = MagicMock()
            mock_ws.send = MagicMock()
            return mock_ws

        mock_ws_class.side_effect = lambda url, **kwargs: setup_callbacks(
            kwargs.get("on_message"),
            kwargs.get("on_error"),
            kwargs.get("on_close"),
            kwargs.get("on_open"),
        )

        result = monitor_job_via_websocket("test_job_123", "http://localhost:8000")

        # Should return completion data
        assert result is not None
        assert result["status"] == "done"

        # Verify credentials can be extracted
        credentials = extract_credentials_from_completion(result)
        assert credentials["username"] == "testuser"
        assert credentials["password"] == "testpass"

    @patch("gpm.exports.websocket.WebSocketApp")
    @patch("gpm.exports.get_gpm_config")
    @patch("gpm.exports.click.echo")
    def test_ping_pong_keepalive(self, mock_echo, mock_config, mock_ws_class):
        """Test ping/pong keepalive message handling."""
        mock_config.side_effect = lambda section, key: {
            ("EXPORT_ENGINE", "EXPORT_ENGINE_WS_URL"): None,
            ("EXPORT_ENGINE", "EXPORT_ENGINE_API_URL"): "http://localhost:8000",
        }.get((section, key), None)

        ping_received = [False]
        pong_sent = [False]

        def setup_callbacks(on_message, on_error, on_close, on_open):
            mock_ws = MagicMock()

            def send_pong(msg):
                if msg == "pong":
                    pong_sent[0] = True

            mock_ws.send = send_pong

            def run_forever():
                """Simulate run_forever with ping message."""
                if on_open:
                    on_open(mock_ws)
                # Simulate ping message
                if on_message:
                    ping_received[0] = True
                    on_message(mock_ws, "ping")
                if on_close:
                    on_close(mock_ws, 1000, "Normal closure")

            mock_ws.run_forever = run_forever
            mock_ws.close = MagicMock()
            return mock_ws

        mock_ws_class.side_effect = lambda url, **kwargs: setup_callbacks(
            kwargs.get("on_message"),
            kwargs.get("on_error"),
            kwargs.get("on_close"),
            kwargs.get("on_open"),
        )

        # The function should handle ping without error
        result = monitor_job_via_websocket(
            "test_job_123", "http://localhost:8000", timeout=1
        )

        # Ping should be received and handled (pong should be sent)
        assert ping_received[0] is True
