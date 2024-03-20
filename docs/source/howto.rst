How to
=======

This how-to tutorial is written according to the normal workflow from demultiplexing, through processing and to the various analyses. Although there are many different options (methods) for each step, the general how-to in GPM is the same. The methods are different in the scripts or files within the methods. Of course you can refer to any step according to your need.


Demultiplex
-----------

If you start with BCL raw data, you should start with:

.. code-block:: shell

   gpm demultiplex --help

Currently, there are the following methods available for demultiplexing:

- bcl2fastq
- cellranger_mkfastq
- cellranger_atac_mkfastq
- evercode_WT

By defining the method for demultiplexing, the template scripts and files will be generated as well as the instructions. 

For example, if you want to run ``cellranger_mkfastq``, you can do:

.. code-block:: shell

   gpm demultiplex --method cellranger_mkfastq \
   --raw /path/to/BCL/folder \
   --output /path/where/new/folder/is/created

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

Before initiating a project, you have to know the followings:

- Project name in the format of **YYMMDD_Name1_Name2_Institute_Application**.
- If available, **PATH/to/FASTQ/**.
- If available, **PATH** to ``project.ini`` from demultiplexing which contains the information of the raw data.
- How you want to process the data (see available methods by ``gpm processing --help``) 

If you have everything ready, you can do:

.. code-block:: shell

   gpm init --from-config /PATH/FASTQ/project.ini \
   --fastq /PATH/FASTQ/FASTQ_FOLDER \
   --name YYMMDD_Name1_Name2_Institute_Application \
   --processing nfcore_RNAseq

This command will:

- Create a new folder with the name, YYMMDD_Name1_Name2_Institute_Application
- Duplicate the previous ``project.ini`` and add new information
- Create a subfolder, **nfcore_RNAseq** and generate the template files and scripts for executing this pipeline

If later you want to add any other processing methods in this project, you can do:

.. code-block:: shell

   gpm processing --fastq /PATH/FASTQ/FASTQ_FOLDER \
   --processing nfcore_miRNAseq \
   /PATH/TO/PROJECT/project.ini

Analysis
--------

After processing the data, now you want to perform some customized analyses according to the experimental design or the initial results. GPM also provides a wide range of analyses ready to use. You can check the help messages by:

.. code-block:: shell

   gpm analysis --help

``--report`` can be generated from our templates according to the application; ``--add`` can specify which analysis method you need and generate the template scripts and files. You can view all the available analyses by:

.. code-block:: shell

   gpm analysis --list project.ini

For example, you have a 3'mRNA-Seq run and want to generate the report and do differential expression analysis, you can do the followings:

.. code-block:: shell

   gpm analysis --report RNAseq \
   --add DGEA_RNAseq \
   project.ini

This command will:

- Create a folder ``analysis``
- Generate a ``Analysis_Report_RNAseq.Rmd`` for rendering a html report
- Create the folder ``analysis/DGEA_RNAseq`` and add the scripts and files needed for this analysis

Then you need to check the files within the analysis folder for learning how to continue the analysis. There might be Rmd or JupyterNotebooks for guiding the analysis.

Export
------

After the analysis is done and now you want to export the data to the clients. GPM provides the command ``export`` for soft-linking everything to the export destination (such as web server) and create the ``.htaccess`` and ``.htpasswd`` for your project.

.. note::  The purpose for soft-linking the files is to avoid duplicating any file or folder.

Please check the help message by:

.. code-block:: shell

   gpm export --help

You should run this command from the root of the project folder where ``project.ini`` is.

.. code-block:: shell

   gpm export --config project.ini \
   --symprefix /mnt/nextgen/
   /mnt/web/var/www/html/data/YYMMDD_Name1_Name2_Institute_Application

This command will do the followings:

- Load all the information in ``project.ini``
- Create the folder in the web server ``/mnt/web/var/www/html/data/YYMMDD_Name1_Name2_Institute_Application``
- Generate ``.htaccess`` in this folder according to your configuration
- Generate a user and its login credential and write into ``.htpasswd``. This user name will be extracted from the folder name ``Name1`` from ``YYMMDD_Name1_Name2_Institute_Application``. However, you can also specify by ``--user``.
- Export the folders according to ``config/export.config`` by symbolic links.

.. note::  ``--symprefix`` is crucial here because it defines how the source files are referred from the export destination to the source. In this example, ``/mnt/nextgen/`` refers to the mounting point of the computational server on the web server.

In case you want to create an empty export folder, you can do the following on the web server where you export your data:

.. code-block:: shell

   gpm export --user myclient YYMMDD_Name1_Name2_Institute_Application

This command will still generate ``.htaccess`` and ``.htpasswd``, but leaves the folder empty for you.

Eventually, you can tar those exporting folders for the users to download.

.. code-block:: shell

   gpm tar-export .

This command needs to be executed in the export project folder (web server) and it will:

- Create ``compressed_tar`` folder
- Iterate through every subfolders except ``compressed_tar`` and compress each subfolder including softlinked files/folders
- The file name of the **tar** files is **Project_name_Subfolder_name**.tar
- ``md5`` file is also generated.

In case you want to re-tar any subfolder, you need to delete that tar file first and redo this step. When you are not sure, you can run the script with ``--dry-run`` to see what is going to happen without actually tarring anything.

Clean
-----

GPM also provides clean command to remove the files or folders which you don't want to archive after the projects closed. The regex patterns of those files/folders are defined in ``config/gpm.ini`` section ``[CLEAN]`` by the key ``PATTERNS``. Please read the help message by:

.. code-block:: shell

   gpm clean --help

You can clean multiple projects at the same time:

.. code-block:: shell

   gpm clean ./2022*

Or you can simulate by ``--dry-run``:

.. code-block:: shell

   gpm clean -v -d ./2022*

Archive
-------

The last stage in the life cycle of a project is archiving, which means to backup the whole project to the archive destination and delete the source files. Please read the help messages:

.. code-block:: shell

   gpm archive --help

You can archive multiple projects at the same time:

.. code-block:: shell

   gpm archive ./2022* /PATH/TO/ARCHIVE/SPACE

If you are not sure how much data will be archived or removed, you can use ``--dry-run``:

.. code-block:: shell

   gpm archive --dry-run --verbose ./2022* /PATH/TO/ARCHIVE/SPACE