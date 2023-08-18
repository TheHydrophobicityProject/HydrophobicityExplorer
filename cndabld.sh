#!/bin/bash
#run within a conda environment with the conda-build package.

# version_number=$(awk -F "\"" '/version/ {print $2}' setup.py)
# md5=$(md5sum dist/*${version_number}.tar.gz | awk '{print $1}')

conda build --output-folder ./conda-out/ ./conda/ -c conda-forge
