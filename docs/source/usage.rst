Usage
=====

.. _installation:

Installation
------------

To use GPM, first install it using pip:

.. code-block:: console

   $ pip install genomicspm

By default, some template data will be copied under `~/gpmdata/`. If you want to change this behavior, please define `GPMDATA`:

.. code-block:: console

   $ export GPMDATA=/preferred/for/gpmdata
   $ pip install genomicspm

You can test your installation by:

.. code-block:: console

   $ gpm --help


Configuration
-------------

GPM will copy the config files and the templates into `~/gpmdata/`. Its structure is as below:

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

Here we would like to introduce you these 5 configuration files under `~/gpmdata/config`. You can find all the explanation as comments in each file.

- `gpm.ini <https://github.com/chaochungkuo/GPM/blob/main/config/gpm.ini>`_: Define the configurations for GPM itself globally.
- `environment.ini <https://github.com/chaochungkuo/GPM/blob/main/config/environment.ini>`_: Define the paths of programs or configuration for the machine where you run GPM.
- `analysis.config <https://github.com/chaochungkuo/GPM/blob/main/config/analysis.config>`_: This file define the group, name, and files for all analyses. This file doesn't need to be modified in most cases, except you want to add new analysis.
- `export.config <https://github.com/chaochungkuo/GPM/blob/main/config/export.config>`_: This file define the behavior for ``gpm export``. This file doesn't need to be modified in most cases, except you want to add new analysis.
- `htaccess <https://github.com/chaochungkuo/GPM/blob/main/config/htaccess>`_: This file will be copied into every export folder by ``gpm export``.

Add user configs
----------------

Any modification on those config files will be overwritten when you install GPM again. In order to keep your changes, you can create a new config file with ".user" at the end of the file name. For example, GPM will load ``gpm.ini.user`` prior ``gpm.ini`` if ``gpm.ini.user`` exists.
