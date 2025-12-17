Usage
=====

.. _installation:

Installation
------------

GPM requires Python 3.12 or higher. We recommend using `UV <https://github.com/astral-sh/uv>`_ for installation:

.. code-block:: console

   $ uv pip install gpm

Or install from source:

.. code-block:: console

   $ git clone https://github.com/chaochungkuo/GPM.git
   $ cd GPM
   $ uv pip install -e .

By default, some template data will be copied under ``~/gpmdata/``. If you want to change this behavior, please define ``GPMDATA``:

.. code-block:: console

   $ export GPMDATA=/preferred/for/gpmdata
   $ uv pip install gpm

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


nextflow.config
---------------

GPM is able to manage both global and pipeline-specific ``nextflow.config`` file. Its approach is the followings:

If "nfcore" is included in the method, GPM will check:

- If no pipeline-specific ``nextflow.config`` exists, GPM will copy one from ``GPMDATA/config/nextflow.config`` to the working processing method.

- If there is a pipeline-specific ``nextflow.config``, GPM will append the content of ``GPMDATA/config/nextflow.config`` to the pipeline-specific ``nextflow.config``.

Note that ``GPMDATA/config/nextflow.config`` can also be customized as ``GPMDATA/config/nextflow.config.user`` according to :ref:`customize_user_configs`. For example, some parameters across the whole server (*max_time* or *max_cpus*) can be defined in ``GPMDATA/config/nextflow.config``. Some pipeline specific parameters can be defined in each method folder in ``nfcore_*/nextflow.config``.