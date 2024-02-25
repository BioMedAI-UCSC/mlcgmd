#!/bin/bash

directory_to_check="$1"

# Check if the directory exists
if [ ! -d "$directory_to_check" ]; then
    echo "Directory $directory_to_check does not exist."
    exit 1
fi

# Find all subdirectories within the given directory
subdirectories=$(find "$directory_to_check" -mindepth 1 -maxdepth 1 -type d)

# Loop through each subdirectory
for subdir in $subdirectories; do
    # Check if the subdirectory is empty
    results_dir="$subdir/result"
    if [ -d "$results_dir" ] && [ -z "$(ls -A "$results_dir")" ]; then
        # Remove the subdirectory
        echo "Removing empty directory with empty 'results' subdirectory: $subdir"
        rm -rf "$subdir"
    fi
done

echo "Done"
