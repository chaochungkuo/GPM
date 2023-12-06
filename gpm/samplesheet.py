import os


def remove_end_slash(path):
    """
    Remove the ending slash in the given path.

    :param path: A path.
    :type path: str
    :return: the modified path
    :rtype: str
    """
    if path.endswith("/"):
        path = path.rstrip("/")
    return path


def generate_samples(STRANDEDNESS, FASTQ_DIR, SAMPLESHEET_FILE,
                     READ1_EXTENSION, READ2_EXTENSION, SINGLE_END,
                     SANITISE_NAME, SANITISE_NAME_DELIMITER,
                     SANITISE_NAME_INDEX):

    strandedness = "unstranded"
    if STRANDEDNESS in ["unstranded", "forward", "reverse"]:
        strandedness = STRANDEDNESS

    fastq_dir_to_samplesheet(
        fastq_dir=FASTQ_DIR,
        samplesheet_file=SAMPLESHEET_FILE,
        strandedness=strandedness,
        read1_extension=READ1_EXTENSION,
        read2_extension=READ2_EXTENSION,
        single_end=SINGLE_END,
        sanitise_name=SANITISE_NAME,
        sanitise_name_delimiter=SANITISE_NAME_DELIMITER,
        sanitise_name_index=SANITISE_NAME_INDEX,
    )


def fastq_dir_to_samplesheet(
    fastq_dir,
    samplesheet_file,
    strandedness="unstranded",
    read1_extension="_R1_001.fastq.gz",
    read2_extension="_R2_001.fastq.gz",
    single_end=False,
    sanitise_name=False,
    sanitise_name_delimiter="_",
    sanitise_name_index=1,
    sc=False,
    r16s=False
):
    def sanitize_sample(path, extension):
        """Retrieve sample id from filename"""
        sample = os.path.basename(path).replace(extension, "")
        if sanitise_name:
            sample = sanitise_name_delimiter.join(
                os.path.basename(path).split(sanitise_name_delimiter)[
                    :sanitise_name_index
                ]
            )
        return sample
