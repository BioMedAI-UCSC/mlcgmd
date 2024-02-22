#!/bin/bash

# Function to delete directory if a subfolder within it is empty
delete_directory_if_subfolder_empty() {
    local directory="$1"
    # Loop through each subdirectory
    for subdir in "$directory"/*; do
        # Check if it's a directory and empty
        if [ -d "$subdir" ] && [ -z "$(ls -A "$subdir")" ]; then
            echo "Deleting $directory because $subdir is empty."
            # Delete the parent directory and exit the function
            rm -rf "$directory"
            return
        fi
    done
    echo "No empty subfolders found in $directory."
}

# Example usage:
delete_directory_if_subfolder_empty "/path/to/parent_directory"
