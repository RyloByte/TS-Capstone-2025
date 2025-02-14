#!/bin/bash

# Check if a file argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <file_with_ids>"
    exit 1
fi

# Input file containing IDs
ID_FILE="$1"

# Check if the file exists
if [ ! -f "$ID_FILE" ]; then
    echo "Error: File '$ID_FILE' not found."
    exit 1
fi

# Read each line from the file and download the corresponding file
while IFS= read -r ID; do
    if [[ -n "$ID" ]]; then
        URL="https://alphafold.ebi.ac.uk/files/AF-${ID}-F1-model_v4.pdb"
        echo "Downloading: $URL"
        wget -q --show-progress -P "$DOWNLOAD_DIR" "$URL"
    fi
done < "$ID_FILE"

echo "Downloads completed. Files are saved in '$DOWNLOAD_DIR'."
