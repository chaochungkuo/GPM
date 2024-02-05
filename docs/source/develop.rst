Development
=====

.. _replaceable_variables:
Replaceable variables in template files
---------------------------------------

You are going to create many files and scripts for ``processing`` and ``analysis``, and many elements in those files can be defined from the configuration files. Making GPM able to handle all these elements is crucial to avoid manual works.

There are two levels of configuration: global and project-wise (see :ref:`two_configuration`). The way to define them are a bit different:

- Global configuration (``gpm.ini`` and ``environment.ini``):
    - Just ues the keys in uppercase. This key in your files will be replaced with the corresponding value by GPM.
    - For example, ``RMD_INSTITUTE_NAME``.

- Project-wise configuration (``project.ini``):
    - Please add **PROJECT_** in front of each key (upper case) defined in ``gpm/gpm.py``.
    - Here is a list of them: **date, name1, name2, institute, application, project.ini, project_path, project_name, project_string, bcl_path, demultiplex_path, fastq_path, fastq_multiqc_path, demultiplex_method, processing_path, processing_method, processing_qc_path, organism, genome_assembly, analysis_path, analysis_types, export_URL, export_user, export_password**
    - For example, **fastq_path** should be **PROJECT_FASTQ_PATH** in your codes.

Please refer to other files in ``analysis`` or ``processing`` to learn how to use them.


How to add a new demultiplexing method?
---------------------------------------

When you want to add a new method for demultiplexing, please follow the steps below:

1. Add a folder under ``demultiplex/`` with a clear and meaningful name. This name must be distinguishable and clear.
2. Add the template scripts and files into this folder with proper variables. Please refer to :ref:`replaceable_variables`.
3. Add this method name into ``GPM_DEMULTIPLEX_METHODS`` in ``gpm.ini``.

Please test it by ``gpm demultiplex --help``.


How to add a new processing method?
-----------------------------------

When you want to add a new method for processing the data, please follow the steps below:

1. Add a folder under ``processing/`` with a clear and meaningful name. This name must be distinguishable and clear.
2. Add the template scripts and files into this folder with proper variables. Please refer to :ref:`replaceable_variables`.
3. Add this method name into ``GPM_PROCESSING_METHODS`` in ``gpm.ini``.

Please test it by ``gpm processing --help``.

How to add a new analysis?
--------------------------

The analysis feature is controlled by two aspects:

- What are the relevant files? Rmd? Jupyter-notebook?
- What is its name and category?

After you can answer the questions above, you can follow the steps below:

1. Add your files under the category folder in ``GPM/analysis/``.
2. Define each file in ``GPM/config/analysis.config``.

Please read through :ref:`replaceable_variables` to learn how to utilize the variables managed by GPM.                                                                                                  

