import os
import glob
import sys


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

    def get_fastqs(extension):
        """
        Needs to be sorted to ensure R1 and R2 are in the same order
        when merging technical replicates. Glob is not guaranteed to produce
        sorted results.
        See also https://stackoverflow.com/questions/6773584/how-is-pythons-glob-glob-ordered
        """
        return sorted(
            glob.glob(os.path.join(fastq_dir, f"*{extension}"),
                      recursive=False)
        )

    read_dict = {}

    # Get read 1 files
    for read1_file in get_fastqs(read1_extension):
        sample = sanitize_sample(read1_file, read1_extension)
        if sample not in read_dict:
            read_dict[sample] = {"R1": [], "R2": []}
        read_dict[sample]["R1"].append(read1_file)

    # Get read 2 files
    if not single_end:
        for read2_file in get_fastqs(read2_extension):
            sample = sanitize_sample(read2_file, read2_extension)
            read_dict[sample]["R2"].append(read2_file)

    # Write to file
    if len(read_dict) > 0:
        out_dir = os.path.dirname(samplesheet_file)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        with open(samplesheet_file, "w") as fout:
            if sc:
                header = ["sample", "fastq_1", "fastq_2"]
            elif r16s: 
                header = ["sampleID", "forwardReads", "reverseReads", "run"]
            else:
                header = ["sample", "fastq_1", "fastq_2", "strandedness"]
            fout.write(",".join(header) + "\n")
            for sample, reads in sorted(read_dict.items()):
                for idx, read_1 in enumerate(reads["R1"]):
                    read_2 = ""
                    if idx < len(reads["R2"]):
                        read_2 = reads["R2"][idx]
                    if sc:
                        sample_info = ",".join([sample, read_1, read_2])
                    elif r16s:
                        sample_info = ",".join([sample, read_1, read_2], "1")
                    else:
                        sample_info = ",".join([sample, read_1, read_2,
                                                strandedness])
                    fout.write(f"{sample_info}\n")
    else:
        error_str = (
            "\nWARNING: No FastQ files found so samplesheet has not been created!\n\n"
        )
        error_str += "Please check the values provided for the:\n"
        error_str += "  - Path to the directory containing the FastQ files\n"
        error_str += "  - '--read1_extension' parameter\n"
        error_str += "  - '--read2_extension' parameter\n"
        print(error_str)
        sys.exit(1)
