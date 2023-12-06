![GitHub License](https://img.shields.io/github/license/chaochungkuo/GPM)
![GitHub all releases](https://img.shields.io/github/downloads/chaochungkuo/GPM/total)
[![Documentation Status](https://readthedocs.org/projects/gpm/badge/?version=latest)](https://gpm.readthedocs.io/en/latest/?badge=latest)

# GPM (Genomic Project Manager)

**GPM (Genomic Project Manager)** is a versatile command-line tool designed for managing and automating bioinformatic workflows. Key features include:

- **End-to-End Workflow**: Seamlessly manage demultiplexing, processing, analysis, interactive reporting, and archiving tasks.
- **Customizable**: Easily adapt GPM to your environment, computer, web service, institute, or specific author requirements.
- **Comprehensive NGS Application Coverage**: Covering a wide range of Next-Generation Sequencing (NGS) applications, including RNAseq, tRNAseq, mRNAseq, ChIPseq, ATACseq, CLIPseq, ampliseq, scRNAseq, scVDJseq, scATACseq, miRNAseq, BWGS, WES, 16S, MAG, and more.
- **Ideal for Bioinformatics Services**: A powerful tool for individuals providing bioinformatic services.

## A Short Demo

After installing and configuring GPM on your computer/server, effortlessly complete your project with the following steps:

### 1. Demultiplexing

Using an RNA-Seq workflow as an example, kick off your NGS project with a BCL folder (BCL_Path):

```bash
gpm demultiplex --raw BCL_Path --output FASTQ_output --method bcl2fastq
```

This command generates essential files for demultiplexing:

- **run_bcl2fastq.sh**: For running bcl2fastq.
- **samplesheet_bcl2fastq.csv**: For adding indices.
- **project.ini**: For storing all project-related information.

### 2. Processing

With your FASTQ files in hand, initiate the processing step:

```bash
gpm init --from-config above/project.ini --fastq path/to/fastqs --name 231231_Chao-Chung_Kuo_UKA_RNAseq
```

This script generates the following files:

```markdown
nfcore_3mRNAseq
├── nextflow.config
└── run_nfcore_3mrnaseq.sh
```

Now, you are ready to run nfcore and process your data.

### 3. Analysis

Generate the analysis report:

```bash
gpm analysis project.ini --application RNAseq
```

An Rmd file, **Analysis_Report_RNAseq.Rmd**, will be added under the _analysis_ folder. Explore available analyses:

```bash
gpm analysis project.ini --list
```

Add a set of files for **DGEA_RNAseq**:

```bash
gpm analysis project.ini --add DGEA_RNAseq
```

With its flexibility and scalability, GPM accelerates the template creation process and encourages code reuse.

For additional tutorials and HowTo guides, refer to the documentation site: [gpm.readthedocs.io](https://gpm.readthedocs.io/)