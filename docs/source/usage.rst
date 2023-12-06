Installation
=====

.. _installation:

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
