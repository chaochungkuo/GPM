import csv
import yaml
import subprocess

# 1. copy DGEA Rmd and modify it according to config.yml

###############################################################################
# Generate the required comparisons to be injected into the markdown files
###############################################################################

path_main_report = "../Analysis_Report_RNAseq.Rmd"
path_main_report_html = "../Analysis_Report_RNAseq.html"
path_DGEA_config = 'DGEA_config.yml'


def inserting_report(text_insert):
    # Inserting lines
    report_scripts = []
    with open(path_main_report) as file:
        for line in file.readlines():
            if "<!-- INSERTING POINT -->" in line:
                report_scripts.append(text_insert)
            report_scripts.append(line)
    with open(path_main_report, "w") as file:
        for line in report_scripts:
            print(line, file=file)


# Load DGEA_config.yml
with open(path_DGEA_config, 'r') as file:
    DGEA_config = yaml.safe_load(file)

# Load contrasts.csv
contrasts_file = "contrasts.csv"
samplesheet_file = "samplesheet.csv"
comparisons = []

with open(contrasts_file, 'r') as file:
    csv_reader = csv.DictReader(file)
    for row in csv_reader:
        comparisons.append(row)

for comparison in comparisons:
    # Iterate each comparison
    com = comparison.split(",")
    # condition_control_treated_blockrep,TRUE,condition,control,treated,batch

    # Define modify dict
    modify_dict = {"TITLEDESCRIPTION": com[0].replace("_", " "),
                   "DGEA_ORGANISM": DGEA_config["default"]["organism"],
                   "DGEA_FILETAG": com[0],
                   "DGEA_VARIABLE": com[2],
                   "DGEA_REF": com[3],
                   "DGEA_TARGET": com[4],
                   "DGEA_spikein_ERCC": DGEA_config["default"]["spikein_ERCC"],
                   "DGEA_PAIRED_BOOLEAN": com[1],
                   "DGEA_countsFromAbundance":
                       DGEA_config["default"]["countsFromAbundance"],
                   "DGEA_lengthcorrection":
                       DGEA_config["default"]["lengthcorrection"]}

    Rmd_template = []
    with open("DGEA_template.Rmd") as file:
        for line in file.readlines():
            for key, value in modify_dict.items():
                if key in line:
                    line = line.replace(key, value)
            if "DGEA_PAIRED" in line:
                if com[1] == "TRUE":
                    line = line.replace("DGEA_PAIRED", "paired")
                elif com[1] == "FALSE":
                    line = line.replace("DGEA_PAIRED", "unpaired")
            Rmd_template.append(line)

    with open("DGEA_"+com[0]+".Rmd", "w") as file:
        for line in Rmd_template:
            print(line, file=file)
    Rmd_name = "DGEA_"+com[0]+".Rmd"
    html_name = "DGEA_"+com[0]+".html"
    commands = ["previous_env=$(conda info --envs | grep '*' | awk '{print $1}')",
                "source activate /opt/miniconda3/envs/rstudio",
                "Rscript -e 'rmarkdown::render('" + Rmd_name + "', output_format = 'html_document', output_file = '"+html_name+"')'",
                "source activate $previous_env"]
    for command in commands:
        subprocess.run(command.split(" "))

inserting_report("# Differential Gene Expression Analysis")
for comparison in comparisons:
    com = comparison.split(",")
    html_name = "DGEA/DGEA_"+com[0]+".html"
    inserting_report("["+com[0].replace("_", " ")+"]("+html_name+")")

commands = ["previous_env=$(conda info --envs | grep "*" | awk '{print $1}')",
            "source activate /opt/miniconda3/envs/rstudio",
            "Rscript -e 'rmarkdown::render('" + path_main_report + "', output_format = 'html_document', output_file = '"+path_main_report_html+"')'",
            "source activate $previous_env"]
for command in commands:
    subprocess.run(command.split(" "))
