How to
=======

This how-to tutorial is written according to the normal workflow from demultiplexing, through processing and to the various analyses. Although there are many different options (methods) for each step, the general how-to in GPM is the same. The methods are different in the scripts or files within the methods. Of course you can refer to any step according to your need.


Demultiplex
-----------

If you start with BCL raw data, you should start with:

.. code-block:: shell

   gpm demultiplex --help

By defining the method for demultiplexing, the template scripts and files will be generated as well as the instructions. For example, if you want to run ``cellranger_mkfastq``, you can do:

.. code-block:: shell

   gpm demultiplex --method cellranger_mkfastq --raw /path/to/BCL/folder --output /path/where/new/folder/is/created

This command will create a folder with the same name as the BCL folder under the defined output folder, and then add the following files in it:
- run_cellranger_mkfastq.sh
- run_merge_lanes.sh
- samplesheet_cellranger.csv

.. note::  The reason to keep the output folder with the same name as the BCL folder is that it is easier for tracing back. In addition, one sequencing run might include the reads from several projects and we have to do the demultiplexing together. This is why we don't want to make demultiplexing step project specific.

These files are everything you need for this task and you need to go through these files and follow the instruction inside to modify it for your need. Then you can run it with:

.. code-block:: shell

   bash run_cellranger_mkfastq.sh

.. note::  It is recommended to run the command within a detachable session such as screen or tmux.

GPM just populates the template files for you. You still need to read and understand how cellranger or bcl2fastq work, and how to define the corresponding samplesheets.

Processing
----------

After you have FASTQ files, you can initiate a new project by ``gpm init``. Please check the help message by:

.. code-block:: shell

   gpm init --help

This command will create the ``project.ini`` and 


Analysis
----------