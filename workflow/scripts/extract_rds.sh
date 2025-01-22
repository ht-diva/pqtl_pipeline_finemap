#!/bin/bash

# Base directory and target directory
BASE_DIR="/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/results/meta_hole_3M/cojo"
TARGET_DIR="/scratch/dariush.ghasemi/projects/rds_conditional_data"

# Create the target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Find and process .rds files
find "$BASE_DIR" -type f -name "*.rds" ! -path "*/finemaping/*" | while read -r file; do

    upstream_dir=$(echo "$file" | grep -oP 'seq\.\d+\.\d+')  # Extract the upstream directory matching seq.*.*
    
    if [[ -n "$upstream_dir" ]]; then
        # Construct the new file name with the upstream directory name
        filename=$(basename "$file")
        new_name="${upstream_dir}_${filename}"
        
        # Copy the file to the target directory
        cp "$file" "$TARGET_DIR/$new_name"
        echo "Copied: $file to $TARGET_DIR/$new_name"
    fi
done


# Compress the final folder into a zip file
zip -r "${TARGET_DIR}.zip" "$TARGET_DIR"
echo "Compressed folder into ${TARGET_DIR}.zip"
