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
gpm init --from-config above/project.ini --fastq path/to/fastqs --name 231231_Chao-Chung_Kuo_UKA_RNAseq --processing nfcore_3mRNAseq
```

This script generates the following files:

```markdown
nfcore_3mRNAseq
├── nextflow.config
└── run_nfcore_3mrnaseq.sh
```

Now, you are ready to run nfcore and process your data.

### 3. Analysis

Generate the analysis report for the given application:

```bash
gpm analysis project.ini --report RNAseq
```

An Rmd file, **Analysis_Report_RNAseq.Rmd**, will be added under the _analysis_ folder. Now you can explore available analyses provided by GPM by:

```bash
gpm analysis project.ini --list
```

Now you can select the analysis you want by name and add a set of files for **DGEA_RNAseq**:

```bash
gpm analysis project.ini --add DGEA_RNAseq
```

Whenever you want to add additional analyses, you can simply repeat the steps above. For example, now we can add **GSEA_hallmarkgenes**.

```bash
gpm analysis project.ini --add GSEA_hallmarkgenes
```

With its flexibility and scalability, GPM accelerates the template creation process and encourages code reuse.

### 4. Export

After you finish the analyses and reports, now you want to export all the results either to a web server or another local folder. No files will be copied, instead, they will be soft linked only.

Below command is used for exporting the results to a new folder in a web server:

```bash
gpm export --config project.ini --symprefix /mnt/nextgen/ --tar /path/to/web/server/folder/
```

Or, you can export the results to another folder locally.

```bash
gpm export --config project.ini --tar /path/to/a/local/folder/
```

### 5. Clean

After you finish the analysis, you might want to clean the temporary files which you don't want to archive. The files for cleaning are defined by regex pattern in `gpm.ini`, section `[CLEAN]` and the key `PATTERNS`,

```bash
gpm clean /path/to/the/folder/
```

### 6. Archive

When everything is done and you want to archive the scripts and codes to another place, you can use the command below:

```bash
gpm archive -v SOURCE_FOLDERS /target/archive/folder/
```

You can use wild cards to define the source folders, such as:

```bash
gpm archive -v 2306* /target/archive/folder/
```


For additional tutorials and HowTo guides, refer to the documentation site: [gpm.readthedocs.io](https://gpm.readthedocs.io/)