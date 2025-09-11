#!/bin/bash
set -euo pipefail

echo "Downloading SA Polony dataset..."
wget -O SA_Polony_100_download.sh \
  https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/SA_Polony_100_download.sh

chmod +x SA_Polony_100_download.sh
bash SA_Polony_100_download.sh

echo "Download finished. Inspect raw/ directory."
