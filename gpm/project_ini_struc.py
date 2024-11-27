from collections import OrderedDict


tags_GPM = OrderedDict([("Project", ["date", "name1", "name2", "institute",
                                     "application", "project.ini",
                                     "project_path", "project_name",
                                     "project_string"]),
                        ("Raw data", ["bcl_path"]),
                        ("Demultiplexing", ["demultiplex_path",
                                            "fastq_path",
                                            "fastq_multiqc_path",
                                            "demultiplex_method"]),
                        ("Processing", ["processing_path",
                                        "processing_method",
                                        "processing_qc_path",
                                        "organism",
                                        "genome_assembly"]),
                        ("Analysis", ["analysis_path",
                                      "analysis_types"]),
                        ("Export", ["export_URL",
                                    "report_URL",
                                    "export_user",
                                    "export_password"]),
                        ("History", ["logs"])])
