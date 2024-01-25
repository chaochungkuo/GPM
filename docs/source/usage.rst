Usage
=====

.. _installation:

Installation
------------

To use GPM, first install it using pip:

.. code-block:: console

   $ pip install genomicspm

By default, some template data will be copied under ``~/gpmdata/``. If you want to change this behavior, please define ``GPMDATA``:

.. code-block:: console

   $ export GPMDATA=/preferred/for/gpmdata
   $ pip install genomicspm

You can test your installation by:

.. code-block:: console

   $ gpm --help


Configuration
-------------

GPM will copy the config files and the templates into ``~/gpmdata/``. Its structure is as below:

.. code-block:: shell

   ~/gpmdata/
   ├── analysis
   ├── config
   │   ├── analysis.config
   │   ├── environment.ini
   │   ├── export.config
   │   ├── gpm.ini
   │   └── htaccess
   ├── demultiplex
   └── processing

Here we would like to introduce you these 5 configuration files under ``~/gpmdata/config``. You can find all the explanation as comments in each file.

- `gpm.ini <https://github.com/chaochungkuo/GPM/blob/main/config/gpm.ini>`_: Define the configurations for GPM itself globally.
- `environment.ini <https://github.com/chaochungkuo/GPM/blob/main/config/environment.ini>`_: Define the paths of programs or configuration for the machine where you run GPM.
- `analysis.config <https://github.com/chaochungkuo/GPM/blob/main/config/analysis.config>`_: This file define the group, name, and files for all analyses. This file doesn't need to be modified in most cases, except you want to add new analysis.
- `export.config <https://github.com/chaochungkuo/GPM/blob/main/config/export.config>`_: This file define the behavior for ``gpm export``. This file doesn't need to be modified in most cases, except you want to add new analysis.
- `htaccess <https://github.com/chaochungkuo/GPM/blob/main/config/htaccess>`_: This file will be copied into every export folder by ``gpm export``.

Workflow
-----------

If you start with BCL raw data, you should start with:

.. code-block:: shell
   gpm demultiplex --help

By defining the method for demultiplexing, the template scripts and files will be generated as well as the instructions.

After you have FASTQ files, you can initiate a new project by ``gpm init``. Please check the help message by:

.. code-block:: shell
   gpm init --help

This command will create the ``project.ini`` and 


nextflow.config
---------------

GPM is able to manage both global and pipeline-specific ``nextflow.config`` file. Its approach is the followings:

If "nfcore" is included in the method, GPM will check:

- If no pipeline-specific ``nextflow.config`` exists, GPM will copy one from ``GPMDATA/config/nextflow.config`` to the working processing method.

- If there is a pipeline-specific ``nextflow.config``, GPM will append the content of ``GPMDATA/config/nextflow.config`` to the pipeline-specific ``nextflow.config``.

Note that ``GPMDATA/config/nextflow.config`` can also be customized as ``GPMDATA/config/nextflow.config.user`` according to :ref:`customize_user_configs`. For example, some parameters across the whole server (*max_time* or *max_cpus*) can be defined in ``GPMDATA/config/nextflow.config``. Some pipeline specific parameters can be defined in each method folder in ``nfcore_*/nextflow.config``.