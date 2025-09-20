#!/bin/bash

# Create directory to store files if it doesn't already exist
mkdir -p clean_files

# List of dataset URLs
urls=(
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/fm_1.fq"
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/fm_2.fq"
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/fm_3.fq"
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/m_4.fq"
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/m_5.fq"
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/m_6.fq"
)

# Download each file if it does not already exist
for url in "${urls[@]}"; do
    wget -nc -P clean_files "$url"
done

echo "✅ All RNA sequencing files downloaded into clean_files/"
