Concept
=====

One project.ini Across All Phases
----------------

In the lifecycle of a project, spanning demultiplexing, processing, analyses, exporting, and archiving, GPM ensures consistent project management by maintaining a single project.ini file across all stages. This file serves as a common thread inherited by subsequent steps, fostering traceability and repeatability for every file and operation.

Keeping Everything in Projects for Repeatability
----------------

Rather than executing functions directly for a project, GPM's approach involves populating the necessary files and scripts specific to defined tasks. Every file within the project guarantees that analyses remain repeatable, even in the absence of GPM.

Flexible Addition of Methods in Any Phase
----------------

GPM is designed to offer flexibility in seamlessly incorporating new methods. Whether it's integrating a demultiplexing method for a new kit, using a custom pipeline for data processing instead of nfcore or cellranger, or adding an R Markdown or Jupyter notebook for future analyses, GPM empowers users to adapt effortlessly to evolving project requirements.