#!/bin/bash

# Directory containing structure files
STRUCTURE_DIR="structures"
OUTPUT_DIR="foldseek_output"
TMP_DIR="temp"

# Ensure the structure directory exists
if [ ! -d "$STRUCTURE_DIR" ]; then
    echo "Error: Directory '$STRUCTURE_DIR' not found."
    exit 1
fi

# Iterate through files in the structure directory
for FILE in "$STRUCTURE_DIR"/*; do
    if [ -f "$FILE" ]; then
        FILENAME=$(basename "$FILE")
        ID=$(echo "$FILENAME" | sed -E 's/AF-([^-]+)-F1-model_v4.pdb/\1/')
        OUTPUT_FILE="fs_easy_search_${ID}.tsv"
        echo "Processing: $FILENAME"
        foldseek easy-search "$FILE" "$STRUCTURE_DIR/" "$OUTPUT_DIR/$OUTPUT_FILE" "$TMP_DIR"
        wait
    fi
done

echo "Processing complete."