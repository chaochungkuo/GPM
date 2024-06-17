
suppressMessages(library(optparse))
suppressMessages(library(clustree))


# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Path for the input file", metavar = "character"),
  make_option(c("-p", "--prefix"), type = "character", default = NULL, help = "the prefix of the clustree", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output file path for the saved dataframe", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$input) || is.null(opt$prefix) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Please provide the paths to the input file, prefix, and the output file.", call. = FALSE)
}


# Load the input file
df <- read.csv(opt$input)
fig = clustree::clustree(df, prefix = "leiden_")
ggsave("clustree.png", fig, width = 10, height = 10, path = opt$output)