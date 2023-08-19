#!/bin/bash
version_number=$(awk -F "\"" '/version/ {print $2}' setup.py)

find ./conda-out/linux-64/ -name "*${version_number}*bz2" -exec conda convert --platform all "{}" -o conda-out \;
find . -name "hydro*ex*${version_number}*bz2" -exec anaconda upload "{}" \;
