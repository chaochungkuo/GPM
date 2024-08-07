#!/bin/bash

# Function to convert a path to absolute path
to_absolute_path() {
    local path="$1"
    if [[ "$path" = /* ]]; then
        echo "$path"
    else
        echo "$(cd "$(dirname "$path")" && pwd)/$(basename "$path")"
    fi
}

# Initialize variables
INPUT_FILE=""
OUTPUT_FILE=""
FROM_FORMAT=""
TO_FORMAT=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -f|--from)
            FROM_FORMAT="$2"
            shift 2
            ;;
        -t|--to)
            TO_FORMAT="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$INPUT_FILE" ] || [ -z "$OUTPUT_FILE" ] || [ -z "$FROM_FORMAT" ] || [ -z "$TO_FORMAT" ]; then
    echo "Usage: $0 -i|--input <input_file> -o|--output <output_file> -f|--from <from_format> -t|--to <to_format>"
    exit 1
fi

# Convert paths to absolute paths
INPUT_FILE=$(to_absolute_path "$INPUT_FILE")
OUTPUT_FILE=$(to_absolute_path "$OUTPUT_FILE")

# Get input and output directories
INPUT_DIR=$(dirname "$INPUT_FILE")
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")

# Get input and output filenames
INPUT_FILENAME=$(basename "$INPUT_FILE")
OUTPUT_FILENAME=$(basename "$OUTPUT_FILE")

# Run Docker container
docker run --rm \
    -v "$INPUT_DIR":/input \
    -v "$OUTPUT_DIR":/output \
    shared-sc_conversion \
    --input "/input/$INPUT_FILENAME" \
    --output "/output/$OUTPUT_FILENAME" \
    --from "$FROM_FORMAT" \
    --to "$TO_FORMAT"

# Check if the command was successful
if [ $? -eq 0 ]; then
    echo "Conversion completed successfully."
    echo "Output file: $OUTPUT_FILE"
else
    echo "Conversion failed."
fi

