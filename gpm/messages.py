import click
from pathlib import Path
from gpm.helper import get_gpm_config

gpm_messages = {
    "demultiplex": {},
    "processing": {},
    "analysis": {}
}

gpm_messages["demultiplex"]["bcl2fastq"] = [
    "1. Modify samplesheet.csv with the proper information. ",
    "   Please add Sample_Project with the correct format ",
    "   (YYMMDD_Name1_Name2_Institute_Application).",
    "2. Check and modify run_bcl2fastq.sh",
    "3. Run run_bcl2fastq.sh with the command below: ",
    "   (Recommend to run it in screen session)",
    "   bash run_bcl2fastq.sh"
]

gpm_messages["demultiplex"]["cellranger_mkfastq"] = [
    "1. Modify samplesheet_cellranger.csv with the proper information. ",
    "2. Check and modify run_cellranger_mkfastq.sh",
    "3. Run run_cellranger_mkfastq.sh with the command below: ",
    "   (Recommend to run it in screen session)",
    "   bash run_cellranger_mkfastq.sh"
]

gpm_messages["demultiplex"]["cellranger_atac_mkfastq"] = [
    "1. Modify samplesheet_cellranger.csv with the proper information. ",
    "2. Check and modify run_cellranger-atac_mkfastq.sh",
    "3. Run run_cellranger-atac_mkfastq.sh with the command below: ",
    "   (Recommend to run it in screen session)",
    "   bash run_cellranger-atac_mkfastq.sh"
]

gpm_messages["processing"]["nfcore_3mRNAseq"] = [
    "1. Generate samplesheet.csv with the following command:",
    "   gpm samplesheet_rnaseq --help",
    "2. Check and modify run_nfcore_3mrnaseq.sh",
    "3. Run run_nfcore_3mrnaseq.sh with the command below: ",
    "   (Recommend to run it in screen session)",
    "   bash run_nfcore_3mrnaseq.sh"
]

gpm_messages["analysis"]["DGEA_RNAseq"] = [
    "1. Generate analysis/samplesheet.csv with all sample information.",
    "   You can try to modify the samplesheet.csv from nfcore pipeline:",
    "   cut -d ',' -f 1 nfcore_RNAseq/samplesheet.csv | \\",
    "   awk 'BEGIN{FS=OFS=\"_\"} {print $0, $1, $2, $3, $4}' OFS=',' \\",
    "   > analysis/samplesheet.csv",
    "2. Modify analysis/DGEA/DGEA_constructor.Rmd and generate sub-reports.",
    "3. Insert the headings with hyperlinks in Analysis_Report.Rmd."
]


class DisplayablePath(object):
    display_filename_prefix_middle = '├──'
    display_filename_prefix_last = '└──'
    display_parent_prefix_middle = '    '
    display_parent_prefix_last = '│   '

    def __init__(self, path, parent_path, is_last):
        self.path = Path(str(path))
        self.parent = parent_path
        self.is_last = is_last
        if self.parent:
            self.depth = self.parent.depth + 1
        else:
            self.depth = 0

    @property
    def displayname(self):
        if self.path.is_dir():
            return self.path.name + '/'
        return self.path.name

    @classmethod
    def make_tree(cls, root, parent=None, is_last=False, criteria=None):
        root = Path(str(root))
        criteria = criteria or cls._default_criteria

        displayable_root = cls(root, parent, is_last)
        yield displayable_root

        children = sorted(list(path
                               for path in root.iterdir()
                               if criteria(path)),
                          key=lambda s: str(s).lower())

        ignore_paths = get_gpm_config("GPM", "GPM_TREE_IGNORE")
        new_children = []
        for child in children:
            # print(child)
            tag_ignore = False
            for p in ignore_paths:
                if p in str(child):
                    tag_ignore = True
            if not tag_ignore:
                new_children.append(child)

        count = 1
        for path in new_children:
            is_last = count == len(new_children)
            if path.is_dir():
                yield from cls.make_tree(path,
                                         parent=displayable_root,
                                         is_last=is_last,
                                         criteria=criteria)
            else:
                yield cls(path, displayable_root, is_last)
            count += 1

    @classmethod
    def _default_criteria(cls, path):
        return True

    def displayable(self):
        if self.parent is None:
            return self.displayname

        _filename_prefix = (self.display_filename_prefix_last
                            if self.is_last
                            else self.display_filename_prefix_middle)

        parts = ['{!s} {!s}'.format(_filename_prefix,
                                    self.displayname)]

        parent = self.parent
        while parent and parent.parent is not None:
            parts.append(self.display_parent_prefix_middle
                         if parent.is_last
                         else self.display_parent_prefix_last)
            parent = parent.parent

        return ''.join(reversed(parts))


def show_tree(target_path):
    click.echo(click.style("\nThe current status in the target directory:",
                           fg='bright_green'))
    paths = DisplayablePath.make_tree(Path(target_path))
    for path in paths:
        click.echo(path.displayable())
    click.echo("")


def show_instructions(command, method):
    click.echo(click.style(" ".join(["Further instructions for",
                                     command,  method+":"]),
                           fg='bright_green'))
    try:
        for line in gpm_messages[command][method]:
            click.echo(line)
    except Exception:
        click.echo("    Not available.")
    click.echo("")
