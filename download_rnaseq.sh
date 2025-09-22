#!/bin/bash

# === Step 1: Create clean_files directory ===
mkdir -p clean_files

# === Step 2: Download RNA-Seq data ===
urls=(
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/fm_1.fq"
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/fm_2.fq"
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/fm_3.fq"
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/m_4.fq"
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/m_5.fq"
  "https://github.com/josoga2/bash-course/raw/refs/heads/main/bash/module_9/START_Here/m_6.fq"
)

for url in "${urls[@]}"; do
    wget -nc -P clean_files "$url"
done

echo "✅ All RNA sequencing files downloaded into clean_files/"

# === Step 3: Git add, commit, and push ===
git add .
git commit -m "Updated RNA-Seq dataset: $(date '+%Y-%m-%d %H:%M:%S')"
git push origin main
