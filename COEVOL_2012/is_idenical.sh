#!/bin/bash

# Get a list of all files in the current directory
files=(nwks_long/*)

# Calculate the MD5 checksum of the first file
first_checksum=$(md5sum "${files[0]}" | awk '{print $1}')

# Iterate through the rest of the files and compare their checksums
for ((i=1; i<${#files[@]}; i++)); do
    checksum=$(md5sum "${files[i]}" | awk '{print $1}')
    if [ "$checksum" != "$first_checksum" ]; then
        echo "Files are not identical: ${files[0]} and ${files[i]}"
        exit 1
    fi
done

echo "All files are identical."
