import os
from pathlib import Path
import pytest
from click.testing import CliRunner
from gpm.main import main
import tempfile
import shutil

def temp_folders():
    tmp_root = tempfile.gettempdir()
    bcl_dir = Path(tmp_root) / "200101_test_0915_TEST/"
    bcl_dir.mkdir(exist_ok=True)
    out_dir = Path(tmp_root) / "output/"
    out_dir.mkdir(exist_ok=True)
    out_demultiplex = os.path.join(out_dir, "200101_test_0915_TEST/")
    return([tmp_root, bcl_dir, out_dir, out_demultiplex])

def rm_output(dir_name):
    if os.path.exists(dir_name) and os.path.isdir(dir_name):
        shutil.rmtree(dir_name)
            
def test_demultiplex_command():
    """Test the demultiplex command."""
    
    runner = CliRunner()
    [tmp_root, bcl_dir, out_dir, out_demultiplex] = temp_folders()
    
    # Test with valid inputs
    rm_output(out_demultiplex)
    result = runner.invoke(main, [
        'demultiplex',
        '--method', 'bcl2fastq',
        '--raw', bcl_dir,
        '--output', out_dir
    ])
    print(result.output)
    assert result.exit_code == 0
    assert "The current status in the target directory" in result.output
    assert "Further instructions for demultiplex bcl2fastq" in result.output
    
def test_demultiplex_invalid_inputs():
    """Test the demultiplex command."""
    
    runner = CliRunner()
    [tmp_root, bcl_dir, out_dir, out_demultiplex] = temp_folders()
    # Test with invalid inputs
    rm_output(out_demultiplex)
    result = runner.invoke(main, [
        'demultiplex',
        '--method', 'bcl2fastq',
        '--raw', "BCL_RAW",
        '--output', out_dir
    ])
    
    assert result.exit_code == 1
    assert "The given path for raw data doesn't exist." in result.output

def test_demultiplex_invalid_methods():
    """Test the demultiplex command."""
    
    runner = CliRunner()
    [tmp_root, bcl_dir, out_dir, out_demultiplex] = temp_folders()
    # Test with invalid method
    rm_output(out_demultiplex)
    result = runner.invoke(main, [
        'demultiplex',
        '--method', 'invalid_method',
        '--raw', out_demultiplex,
        '--output', out_dir
    ])
    
    assert result.exit_code != 1
    assert "Invalid value" in result.output

def test_init_project_command():
    """Test project initialization command."""
    runner = CliRunner()
    [tmp_root, bcl_dir, out_dir, out_demultiplex] = temp_folders()
    # Test basic project initialization
    # with runner.isolated_filesystem(temp_dir=tmp_root):
    old_cwd = os.getcwd()
    os.chdir(tmp_root)
    project_dir = Path(tmp_root) / "240101_Test_Test_Institute_RNAseq"
    rm_output(project_dir)
    
    result = runner.invoke(main, [
        'init',
        '--name', '240101_Test_Test_Institute_RNAseq'
    ])
    
    assert result.exit_code == 0
    print(result.output)
    assert project_dir.exists()
    assert (project_dir / "project.ini").exists()
    os.chdir(old_cwd)
    
# def test_processing_command(tmp_project_dir, mock_project_config, mock_fastq_data):
#     """Test the processing command."""
#     runner = CliRunner()
    
#     # First initialize a project
#     runner.invoke(main, [
#         'init',
#         '--name', 'TestProject',
#         '--fastq', str(mock_fastq_data)
#     ])
    
#     # Test processing command
#     result = runner.invoke(main, [
#         'processing',
#         str(mock_project_config),
#         '--fastq', str(mock_fastq_data),
#         '--processing', 'RNAseq'
#     ])
    
#     assert result.exit_code == 0
#     assert "Processing completed" in result.output

# def test_analysis_command(tmp_project_dir, mock_project_config, mock_analysis_template):
#     """Test the analysis command."""
#     runner = CliRunner()
    
#     # Test listing available templates
#     result = runner.invoke(main, [
#         'analysis',
#         str(mock_project_config),
#         '--list'
#     ])
    
#     assert result.exit_code == 0
#     assert "Available templates" in result.output
    
#     # Test adding a template
#     result = runner.invoke(main, [
#         'analysis',
#         str(mock_project_config),
#         '--add', 'template1'
#     ])
    
#     assert result.exit_code == 0
#     assert "Template added" in result.output

# def test_report_command(tmp_project_dir, mock_project_config, mock_report_template):
#     """Test the report generation command."""
#     runner = CliRunner()
    
#     result = runner.invoke(main, [
#         'report',
#         str(mock_project_config),
#         '--report', 'RNAseq'
#     ])
    
#     assert result.exit_code == 0
#     assert "Report generated" in result.output
    
#     # Verify report was created
#     report_dir = tmp_project_dir / "reports"
#     assert (report_dir / "report.html").exists()

# def test_samplesheet_commands(tmp_project_dir, mock_fastq_data):
#     """Test samplesheet generation commands."""
#     runner = CliRunner()
    
#     # Test RNA-seq samplesheet generation
#     samplesheet_path = tmp_project_dir / "samplesheet.csv"
#     result = runner.invoke(main, [
#         'samplesheet_rnaseq',
#         str(samplesheet_path),
#         str(mock_fastq_data),
#         '--st', 'unstranded'
#     ])
    
#     assert result.exit_code == 0
#     assert samplesheet_path.exists()
    
#     # Test scRNA-seq samplesheet generation
#     sc_samplesheet_path = tmp_project_dir / "sc_samplesheet.csv"
#     result = runner.invoke(main, [
#         'samplesheet_scrnaseq',
#         str(sc_samplesheet_path),
#         str(mock_fastq_data)
#     ])
    
#     assert result.exit_code == 0
#     assert sc_samplesheet_path.exists()

# def test_export_command(tmp_project_dir, mock_project_config):
#     """Test the export command."""
#     runner = CliRunner()
    
#     export_dir = tmp_project_dir / "export"
#     result = runner.invoke(main, [
#         'export',
#         str(export_dir),
#         '--config', str(mock_project_config)
#     ])
    
#     assert result.exit_code == 0
#     assert "Export completed" in result.output
#     assert export_dir.exists()

# def test_clean_command(tmp_project_dir):
#     """Test the clean command."""
#     runner = CliRunner()
    
#     # Create some test directories to clean
#     test_dirs = [
#         tmp_project_dir / "dir1",
#         tmp_project_dir / "dir2"
#     ]
#     for d in test_dirs:
#         d.mkdir()
#         (d / "test.txt").touch()
    
#     result = runner.invoke(main, [
#         'clean',
#         str(test_dirs[0]),
#         str(test_dirs[1])
#     ])
    
#     assert result.exit_code == 0
#     assert "Cleaning completed" in result.output

# def test_archive_command(tmp_project_dir):
#     """Test the archive command."""
#     runner = CliRunner()
    
#     # Create source directories
#     source_dirs = [
#         tmp_project_dir / "source1",
#         tmp_project_dir / "source2"
#     ]
#     for d in source_dirs:
#         d.mkdir()
#         (d / "test.txt").touch()
    
#     archive_dir = tmp_project_dir / "archive"
#     result = runner.invoke(main, [
#         'archive',
#         str(source_dirs[0]),
#         str(source_dirs[1]),
#         str(archive_dir)
#     ])
    
#     assert result.exit_code == 0
#     assert "Archiving completed" in result.output
#     assert archive_dir.exists() 