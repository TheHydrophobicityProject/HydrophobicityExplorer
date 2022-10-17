#!/bin/bash
#run within a conda environment with the conda-build package.
conda build --output-folder ./conda-out/ ./conda/ -c conda-forge
