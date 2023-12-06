Installation
=====

.. _installation:

To use GPM, first install it using pip:

.. code-block:: console

   (.venv) $ pip install genomicspm

By default, some template data will be copied under ``~/gpmdata/``. If you want to change this behavior, please define ``GPMDATA``:

.. code-block:: console

   (.venv) $ export GPMDATA=/preferred/for/gpmdata
   (.venv) $ pip install genomicspm

You can test your installation by:

.. code-block:: console

   (.venv) $ gpm --help
