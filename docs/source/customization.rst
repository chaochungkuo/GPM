Customization
=====

.. _customize_user_configs:
Customize your own config files
----------------
Any modification on those config files will be overwritten when you install GPM again. In order to keep your changes, you can create a new config file with ".user" at the end of the file name. For example, GPM will load ``gpm.ini.user`` prior ``gpm.ini`` if ``gpm.ini.user`` exists.

Customize your own author information
----------------
You can edit ``gpm.ini`` or ``gpm.ini.user`` on the lines below:

.. code-block:: shell

    [RMD]
    RMD_AUTHORS = "First name Last name, Institute, email@address"
                  "First name Last name, Institute, email@address"

You can add multiple persons with proper indentation.

Customize your own institute information
----------------
You can edit ``gpm.ini`` or ``gpm.ini.user`` on the lines below:

.. code-block:: shell

    [RMD]
    RMD_INSTITUTE_NAME = Genomic Facility in IZKF, RWTH Aachen Uniklinik
    RMD_INSTITUTE_LOGO = /path/to/RWTH_IZKF_GF_Logo_rgb.png

