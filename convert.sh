#!/bin/bash
read -rp "Enter path to the system-specific package: " path
conda convert --platform all "$path" -o conda-out

version=$(grep version setup.py | cut -d "\"" -f2)
find . -name "hydro*ex*${version}*bz2" -exec anaconda upload "{}" \;
