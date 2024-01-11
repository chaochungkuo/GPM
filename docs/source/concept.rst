Concept
=====

One project.ini across all phases
----------------

In the lifecycle of a project, spanning demultiplexing, processing, analyses, exporting, and archiving, GPM ensures consistent project management by maintaining a single project.ini file across all stages. This file serves as a common thread inherited by subsequent steps, fostering traceability and repeatability for every file and operation.

Two levels of configuration files
----------------

There are two levels of configuration files in GPM:
- Global configuration under `GPMDATA/config/`
    - `environment.ini` by default or `environment.ini.user` if it exists.
    - `gpm.ini` by default or `gpm.ini.user` if it exists.
- Project-wise configuration as `project.ini`

All the options defined in the above files are used as the key for inserting their corresponding values while copying any files in GPM in creating reports or adding analysis templates. This offers a flexible way to customize your own script and reports.

Keeping everything in projects for repeatability
----------------

Rather than executing functions directly for a project, GPM's approach involves populating the necessary files and scripts specific to defined tasks. Every file within the project guarantees that analyses remain repeatable, even in the absence of GPM.

Flexible addition of methods in any phase
----------------

GPM is designed to offer flexibility in seamlessly incorporating new methods. Whether it's integrating a demultiplexing method for a new kit, using a custom pipeline for data processing instead of nfcore or cellranger, or adding an R Markdown or Jupyter notebook for future analyses, GPM empowers users to adapt effortlessly to evolving project requirements.